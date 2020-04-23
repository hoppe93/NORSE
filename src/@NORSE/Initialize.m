function Initialize(o, varargin)
    % Initializes variables and sets up parameters and flags.
    % Usage: 
    %   Initialize()
    %   Initialize(useExternalInput)
    %
    % For case 4 of the initialDistribution parameter an externally given
    % distribution and two grid vectors must be provided in a
    % structure form, with the the three inputs being three separate fields
    % of the named structure. For more detail please check case 4 of the
    % initialDistribution initialization.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Print some info
    if ~o.silent    
        fprintf('\n*********************************************************************\n');
        fprintf('                    Starting NORSE calculation\n');
        fprintf('*********************************************************************\n');
        fprintf(' Numerical parameters:\n');
        fprintf('       nP: %d, nXi: %d, nL: %d, pMax: %.2f, dt: %.3g, tMax: %.3g.\n',...
                                            o.grid.nP,o.grid.nXi,o.nL,o.grid.pMax,o.dt,o.tMax);
        fprintf(' Settings:\n');
        if o.pGridMode
            pType = 'nonuniform';
        else
            pType = 'uniform';
        end
        if o.xiGridMode
            xiType = 'nonuniform';
        else
            xiType = 'uniform';
        end
        fprintf('   Using a %s grid in p and a %s grid in xi.\n',pType,xiType);   
        switch o.runawayRegionMode
            case 0
                reType = 'isotropic force balance';
            case 1
                reType = 'xi-dependent force balance';
            case 2
                reType = 'linear momentum-space trajectories';
            case 3
                reType = 'non-linear momentum-space trajectories';
            otherwise
                error('Invalid runaway region mode.');
        end
        fprintf('   Calculating runaway region using %s.\n',reType);
    end

    tStartInit = tic;

    
    % Initialize a struct to time different parts of the calculation
    o.timing = struct('initialization',0,...
                      'potentialMats',0,...
                      'timeAdvance',0,...
                      'collOp',0,...
                      'matrixInversion',0,...
                      'matrixFactorization',0,...
                      'GMRES',0,...
                      'postProcess',0,...
                      'separatrix',0,...
                      'total',0);
                  
    
    %%% Handle the physical parameters
    o.ProcessPhysicalParameters();
            
    
    %%% Calculate Legendre polynomials and a matrix for mapping 
    %%% from the grid to the Legendre polynomial representation            
    o.BuildMappingMatrices();            


    %%% Initialize operator matrices to calculate the potenials 
    %%% as functions of f
    o.potentials = Potentials(o);
    tStartPotMat = tic;
    o.potentials.GeneratePotentialMatrices(); 
    o.timing.potentialMats = toc(tStartPotMat);            

    
    %%% Time-advancement and collision operator objects    
    o.InitializeTAandC();
    o.collisionOperator.InitializeParts();    


    %%% Assemble the other components of the matrix               
    o.AssembleEFieldOperator();                                                
    o.AssembleSynchrotronOperator();

    %Initialize the particle source before the heat sink!
    o.particleSource = ParticleSource(o);
    o.heatSink = HeatSink(o);            

    %Initialize save array for source magnitudes       
    s = o.heatSink.magnitudeStructTemplate;
    o.energyChangeMagnitudes = s;
    o.energyChangeMagnitudes(o.nSaveSteps) = s;
    o.densityChangeMagnitudes(o.nSaveSteps) = 0;

    %We will also need an identity matrix for the time advance.
    %Here we can also put boundary conditions.
    I = speye(o.grid.matSize);
    I(sub2ind(size(I),o.grid.idsPMax,o.grid.idsPMax)) = 0; 
                            %Make sure we do not interfere with the
                            %boundary condition at p=p_max
    I(end,end) = 0; %Point at p=0 (where we want ddp = 0)             
    o.identityWithBoundaryConditions = I;            


    %%% Initialize a struct to keep track of the runaway region
    mask = zeros(o.grid.matSize,1);
    pcs = Inf*ones(o.grid.nXi,1);
    o.runawayRegion = struct('mask',mask,'pcs',pcs,'EHat',0,'time',0);
    %Check that settings are consistent
    if o.abortWhenOnlyRunaways && o.runawayRegionMode~=3
        warning('abortWhenOnlyRunaways requires runawayRegionMode=3 and has been disabled.');
        o.abortWhenOnlyRunaways = 0;
    end
    if o.runawayRegionMode==3 && ~o.savePotentials
	warning('runawayRegionMode=3 requires that the potentials be saved. savePotentials has been enabled.')
	o.savePotentials = 1;
    end

    %%% Initialize a distribution 
    %Allocate arrays
    o.f = zeros(o.grid.matSize,o.nSaveSteps);
    if o.savePotentials
        o.Upsilon0 = zeros(o.grid.matSize,o.nSaveSteps);    
        o.Upsilon1 = zeros(o.grid.matSize,o.nSaveSteps);    
        o.Upsilon2 = zeros(o.grid.matSize,o.nSaveSteps); 
        o.Pi = zeros(o.grid.matSize,o.nSaveSteps);
    end

    %Set the initial state
    o.Print('   The starting distribution is ');
    switch o.initialDistribution
        case 0 
            %A relativisitc Maxwellian (in Legendre space)
            f0 = o.maxwellianPreFactor(0) * ...
                    exp( (1-o.grid.gamma)/o.Theta(0) );

            %f(p=0) is treated separately
            initialF = kron([1;zeros(o.nL-1,1)],f0(2:end)); %All other modes=0                        
            initialF = [initialF;f0(1)];                    
            o.f(:,1) = o.MapLegModesToBigVector(initialF); %Map it onto the grid
            o.Print('a Maxwellian.\n');
        case 1
            %A weighted Maxwellian
            delta2 = 2*o.Theta(0);
            f0 = o.grid.p.^2/delta2.*o.maxwellianPreFactor(0) * ...
                                exp( (1-o.grid.gamma)/o.Theta(0) ); 

            initialF = kron([1;zeros(o.nL-1,1)],f0(2:end));                         
            initialF = [initialF;f0(1)];                     
            o.f(:,1) = o.MapLegModesToBigVector(initialF);
            o.Print('a weighted Maxwellian.\n');
        case 2  
            %A shifted Maxwellian   
            pShift = o.maxwellianShift*sqrt(2*o.Theta(0));
            gamma0 = sqrt(1+pShift^2);
            gamma = sqrt(1+o.grid.p2D.^2);
            pPara = o.grid.pPara2D;
            f0 = exp( (1-gamma0*gamma+pShift*pPara)/o.Theta(0) );                    
            o.f(:,1) = o.grid.MapGridToBigVector(...
                                    o.maxwellianPreFactor(0)*f0);
            o.Print('a shifted Maxwellian.\n');
        case 3  
            %Two shifted Maxwellians                    
            yShift = o.maxwellianShift;

            pShift1 = yShift*sqrt(2*o.Theta(0));
            pShift2 = -yShift*sqrt(2*o.Theta(0));
            gamma0 = sqrt(1+pShift1^2);
            gamma = sqrt(1+o.grid.p2D.^2);
            pPara = o.grid.pPara2D;
            f0 = exp( (1-gamma0*gamma+pShift1*pPara)/o.Theta(0) )...
               + exp( (1-gamma0*gamma+pShift2*pPara)/o.Theta(0) );
            o.f(:,1) = o.grid.MapGridToBigVector(...
                                    o.maxwellianPreFactor(0)*f0);
            o.Print('two shifted Maxwellians.\n');
        case 4
            % External distribution
            % As an input, a structure must be given, with three fields
            % corresponding to f, extPBig, extXiBig.
            
            % Check if the external grid is the same as the created
            % Set a limit to the difference
            limit = 1e-12;

            % Check if the gird dimensions ar the same
            pDim = abs(varargin{1}.g.nP - o.grid.nP);
            xiDim = abs(varargin{1}.g.nXi - o.grid.nXi);
            maxDim = abs(varargin{1}.g.pMax - o.grid.pMax);

            dimensions = max([pDim, xiDim, maxDim]);

            if dimensions == 0
                % If the dimesnions are the same, check if the grid points are the identical
                identicalGrid = max(abs((varargin{1}.extPBig-o.grid.pBig)./varargin{1}.extPBig))<limit...
                   &&  max(abs((varargin{1}.extXiBig-o.grid.xiBig)./varargin{1}.extXiBig))< limit;
            else
                identicalGrid = 0;
            end
            
            if identicalGrid
                % If grid is identical, use the external dstribution directly
                if nargin == 2
                    o.f(:,1) = varargin{1}.f;
                    o.Print('an externally given distribution.\n');
                else
                    error('Invalid number of input arguments for external distribution.');
                end
            else
                % Interpolate distribution to the initialized NORSE grid otherwise
                f2D = interp2(varargin{1}.g.xi2D, varargin{1}.g.p2D, varargin{1}.f2D, o.grid.xi2D, o.grid.p2D);
                f = o.grid.MapGridToBigVector(f2D);
                o.f(:,1) = f;
                o.Print('an externally given distribution interpolated to the NORSE grid.\n');
            end
        otherwise
            error('Invalid initialDistribution parameter');
    end        

    %Calculate the potentials from f and save the the initial state
    f = o.f(:,1);
    fls  = o.MapBigVectorToLegModes(f);
    o.potentials.Update(fls);
    o.collisionOperator.Assemble(0);
    o.heatSink.CalculateEnergyChangeMagnitude(f,0,0);                
    o.timeAdvance.SaveStepData(1,f,0);

    
    o.timing.initialization = toc(tStartInit);
end
