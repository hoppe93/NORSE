function InitializeGrid(o)
    % Function that constructs the grid using the parameters and
    % settings set as object properties, and calculates all the
    % associated quantities.
    %
    % Usage:
    %   InitializeGrid()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nP = o.nP;
    nXi = o.nXi;
    o.matSize = (nP-1)*nXi+1;
    o.idsPMax = (nP-1)*(1:nXi);

    %%% Initialize grid in p %%%
    if o.pGridMode
        %A nonuniform grid. 

        %Get a uniform grid on [0,1].
        [s,sWeights,dds,d2ds2] = o.CreateUniformGrid(nP,0,1);

        %Define the nonuniform mapping. Additional grid mappings
        %can be implemented here as long as they define the
        %quantities p, dpds and d2pds2.
        switch o.pGridMode                 
            case 1 %A grid with p = s^2 + pGridParameter*s
                p = s.*s+o.pGridParameter*s;
                dpds = 2*s+o.pGridParameter;
                d2pds2 = 2*ones(size(s));                        
            case 2 %A grid with p = s^3 + pGridParameter*s
                p = s.^3+o.pGridParameter*s;
                dpds = 3*s.^2+o.pGridParameter;
                d2pds2 = 6*s;                        
            case 3 %A grid with p = s^4 + pGridParameter*s
                p = s.^4+o.pGridParameter*s;
                dpds = 4*s.^3+o.pGridParameter;
                d2pds2 = 12*s.^2;                        
            case 4 %A grid with a step in spacing, giving a dense 
                   %grid at low p and a sparse grid at high p                        
                [p,dpds,d2pds2] = o.DefineStepMapping(s);                              
            otherwise
                error('Invalid mode for the p grid.')
        end

        %Remap derivatives and quadrature weights
        ddp = spdiags(1./dpds,0,nP,nP)*dds;
        d2dp2 = -spdiags(d2pds2 ./ (dpds.^3),0,nP,nP) * dds ...
              + spdiags((1./dpds).^2,0,nP,nP)*d2ds2;
        pWeights = dpds .* sWeights;

        %Rescale to get the desired pMax
        scaleFactor = o.pMax/p(end);
        o.p = scaleFactor*p;
        o.ddp = ddp / scaleFactor;
        o.d2dp2 = d2dp2 / (scaleFactor*scaleFactor);
        o.pWeights = pWeights*scaleFactor;
        o.dp = p(2)-p(1);                
    else
        %A uniform grid
        [o.p,o.pWeights,o.ddp,o.d2dp2] = ...
                                o.CreateUniformGrid(nP,0,o.pMax);
        o.dp = o.p(2)-o.p(1);                      
    end        

    %%% Initialize grid in xi %%%
    if o.xiGridMode
        %A nonuniform grid.
        %Additional grid mappings can be implemented here as long 
        %as they define the quantities needed for the remapping below.
        switch o.xiGridMode
            case 1 
                %Use a sigmoid -- denser grid close to xi=1 (and to
                %a lesser extent xi=-1). Specify the grid shape
                %wanted, then iterate to find the exact parameters
                %that give a grid point at xi=0, which we need for
                %the boundary conditions on p.

                c=2.0011; d=16; h=1.001; 
                                %c, which should be >h+1, controls 
                                %the amount of nonlinearity
                [c,h] = o.FindOptimalSigmoidMapping(c,d,h,nXi);

                %Create a uniform grid
                lowerSBound = log((h-1)/(c-h+1))/d;
                upperSBound = log((h+1)/(c-h-1))/d;
                [s,sWeights,dds,d2ds2] = o.CreateUniformGrid(nXi,...
                                            lowerSBound,upperSBound);

                %Define the mapping                        
                o.xi = c./(1+exp(-d*s))-h;

                %Find (and fix) the 0 point
                [~,o.xi0Id] = min(abs(o.xi));
                o.xi(o.xi0Id) = 0;

                %Derivatives of the mapping
                dxids = c*d*exp(d*s)./(1+exp(d*s)).^2; 
                d2xids2 = c*d*d*exp(d*s).*(1-exp(d*s))./(1+exp(d*s)).^3;        
            case 2 
                %Uniform grid in theta = arccos(xi). Improves the
                %resolution close to xi=\pm 1. The mapping cannot
                %be exactly uniform to avoid a vanishing derivative
                %at the endpoints.

                c = 963/610; %Almost pi/2, but not quite, to avoid div by 0
                almostOne = 0.995;
                almostZero = 1e-3;

                %Get a unform grid
                [s,sWeights,dds,d2ds2] = o.CreateUniformGrid(nXi,...
                                              -almostOne,almostOne);

                %Define the mapping and fix the 0 point
                o.xi = sin(c*s);
                o.xi0Id = find(o.xi==0,1);                
                o.xi(o.xi0Id) = 0;

                %Derivatives
                dxids = c*cos(c*s);
                dxids(1) = almostZero;
                dxids(end) = almostZero;
                d2xids2 = -c^2*sin(c*s); 
                d2xids2(o.xi0Id) = almostZero;
            otherwise
                error('Invalid mode for the xi grid.')
        end

        %Remap
        o.ddxi = spdiags(1./dxids,0,nXi,nXi)*dds;
        o.d2dxi2 = -spdiags(d2xids2 ./ (dxids.^3),0,nXi,nXi) ...
                 * dds + spdiags((1./dxids).^2,0,nXi,nXi)*d2ds2;        
        o.xiWeights = dxids .* sWeights;                             
    else
        %A uniform grid
        [o.xi,o.xiWeights,o.ddxi,o.d2dxi2] = o.CreateUniformGrid(...
                                                          nXi,-1,1);      
        %Find the 0 point.
        o.xi0Id = find(o.xi==0,1);                
    end            
    if isempty(o.xi0Id)
        error('The xi grid must have a point at \xi=0. Choose an odd value for nXi.');
    end                

    %%% Calculate aggregate quantities
    o.gamma = sqrt(1+o.p.*o.p);
    [o.xi2D,o.p2D] = meshgrid(o.xi,o.p); % xi_i varies on rows, 
                                         % p_i varies in columns
    o.pBig = o.MapGridToBigVector(o.p2D); 
    o.xiBig = o.MapGridToBigVector(o.xi2D);
    o.gammaBig = sqrt(1+o.pBig.*o.pBig);
    o.pPara2D = o.p2D.*o.xi2D;
    o.pPerp2D = o.p2D.*sqrt(1-o.xi2D.^2); 

    o.ConstructBigDifferentiationMatrices();
    end
