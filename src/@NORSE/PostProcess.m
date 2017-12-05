function PostProcess(o,varargin)
    % Calculates moments of the distribution and other
    % post-processing work. 
    % Usage:
    %   PostProcess()
    %   PostProcess(step)   
    % 
    % If an argument is passed, it specifies the first time step to
    % process (useful for restarts).
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tStartProcess = tic;

    c = o.constants.c;
    e = o.constants.e;                      
    currentNormFactor = e*c*o.fM0;
    energyNormFactor = o.constants.m*c^2*o.fM0;

    if nargin >= 2 && isnumeric(varargin{1}) && isscalar(varargin{1})
        %This is a restart -- we need to extend arrays
        firstId = varargin{1};
        nSteps = o.nSaveSteps-firstId+1; 
        tSep = o.timing.separatrix;
    else
        %This is a new run -- we need to initialize arrays
        firstId = 1;
        nSteps = o.nSaveSteps;
        tSep = 0;
    end

    %%% Calculate moments of f
    %Initialize arrays
    density                  = zeros(1,nSteps); 
    energy                   = zeros(1,nSteps);
    momentum                 = zeros(1,nSteps);
    currentDensity           = zeros(1,nSteps);            
    effectiveTemp            = zeros(1,nSteps);
    effectiveBulkTemperature = zeros(1,nSteps);
    bulkDensity              = zeros(1,nSteps);
    runawayFraction          = zeros(1,nSteps);
    runawayEnergyFraction    = zeros(1,nSteps);

    %Loop throught the time steps and calculate moments
    for i = firstId:o.nSaveSteps
        iSave = i-firstId+1;                
        density(iSave)  = o.grid.intdpdxi'*o.f(:,i);                
        energy(iSave)   = o.grid.intdpdxi'*((o.grid.gammaBig-1).*o.f(:,i));
        momentum(iSave) = o.grid.intdpdxi'...
            *(o.grid.xiBig.*o.grid.pBig.*o.f(:,i));
        currentDensity(iSave)  = o.grid.intdpdxi'...
            *(o.grid.xiBig.*o.grid.pBig./o.grid.gammaBig.*o.f(:,i));       
        effectiveTemp(iSave) = o.CalculateEffectiveTemperature(...
                                        energy(iSave),density(iSave));

        tSepStart = tic;
        [runawayFraction(iSave),runawayEnergyFraction(iSave)]... 
                                    = o.CalculateRunawayFraction(i);
        tSep = tSep + toc(tSepStart);    

        bulkDensity(iSave) = (1-runawayFraction(iSave))*density(iSave);
        bulkEnergy = (1-runawayEnergyFraction(iSave))*energy(iSave);
        effectiveBulkTemperature(iSave) = o.CalculateEffectiveTemperature(...
                                    bulkEnergy,bulkDensity(iSave));
    end      

    %Save the arrays
    if nargin >= 2  
        %Restart -- we need to extend the existing arrays
        o.density  = [o.density,o.fM0*density];
        o.energy   = [o.energy,energyNormFactor*energy];
        o.energyPerParticle = [o.energyPerParticle,...
            energyNormFactor*energy./(o.constants.e*o.fM0*density)];                    
        enM = 3*e/2 * o.n*o.T;
        o.energyOfNonRelMaxwellian = ...
            [o.energyOfNonRelMaxwellian,enM(o.times(firstId:end))];
        o.effectiveTemperature = [o.effectiveTemperature,effectiveTemp];
        o.effectiveBulkTemperature = ...
                [o.effectiveBulkTemperature,effectiveBulkTemperature];
        o.bulkDensity = [o.bulkDensity,o.fM0*bulkDensity];
        o.momentum = [o.momentum,o.fM0*momentum];
        o.currentDensity  = [o.currentDensity,currentNormFactor*currentDensity];
        o.runawayFraction = [o.runawayFraction,runawayFraction];
        o.runawayEnergyFraction = ...
                    [o.runawayEnergyFraction,runawayEnergyFraction];
        o.runawayGrowthRate = o.CalculateRunawayGrowthRate();
        o.isRestart = [o.isRestart,true,false(1,nSteps-1)];
        if o.timeAdvanceMode ~= 3 
            o.EHatInductive = [o.EHatInductive,o.EHat(o.times(firstId:end))];
        end
    else
        %New run                
        o.density  = o.fM0*density;
        o.energy   = energyNormFactor*energy;
        o.energyPerParticle = o.energy./(o.constants.e*o.density);
        enM = 3*e/2 * o.n*o.T;
        o.energyOfNonRelMaxwellian = enM(o.times);
        o.effectiveTemperature = effectiveTemp;
        o.effectiveBulkTemperature = effectiveBulkTemperature;
        o.bulkDensity = o.fM0*bulkDensity;
        o.momentum = o.fM0*momentum;
        o.currentDensity  = currentNormFactor*currentDensity;
        o.runawayFraction = runawayFraction;
        o.runawayEnergyFraction = runawayEnergyFraction;
        o.runawayGrowthRate = o.CalculateRunawayGrowthRate();
        o.isRestart = false(1,nSteps);
        if o.timeAdvanceMode ~= 3 
            o.EHatInductive = o.EHat(o.times);
        end

        o.plot = Plot(o); %Generate a Plot object that we can later use for visualization
    end

    o.timing.postProcess = o.timing.postProcess + toc(tStartProcess);
    o.timing.separatrix  = tSep;
end
