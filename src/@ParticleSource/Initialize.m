function Initialize(o)
    % Initializes the quantities necessary for the particle
    % source/sink.
    %
    % Usage: 
    %   Initialize()            
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~o.norse.n.isScalar 
        %Density changes will have to be made -- make
        %time-dependent parameters for aP and the source density
        %moment
        Theta = o.norse.Theta;
        o.energyContentWBar = besselk(3,1/Theta,1)/besselk(2,1/Theta,1) ...
                                                - 1 - Theta;
        o.particleSourceAP = 2/Theta - 3 ...
                             - (3+3*Theta)/o.energyContentWBar;
        o.particleSourceDensityMoment = ...
                        o.norse.n*( o.energyContentWBar/Theta ...
                                            + o.particleSourceAP );
        o.gammaMin1 = o.grid.gammaBig-1;
    end
    %Initialize an empty operator. This will be overwritten later
    %if needed (otherwise we only need to generate it once)
    o.particleSourceOperator = sparse(o.grid.matSize,1);
end
