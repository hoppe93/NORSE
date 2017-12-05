function SynchSettings(o)
    % Makes sure settings are consistent with those in NORSE.
    % Should be called before the calculation starts (or at restart).
    %
    % Usage:
    %   SynchSettings()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    oN = o.norse;
    
    %Prepare the heat and particle sources
    oN.particleSource.Initialize();
    oN.heatSink.SynchSettings();
    oN.collisionOperator.SynchSettings();
    
    %Synch time-advance settings from the NORSE object
    o.dt                          = oN.dt;
    o.tMax                        = oN.tMax;
    o.nTimeSteps                  = oN.nTimeSteps;
    o.nSaveSteps                  = oN.nSaveSteps;
    o.times                       = oN.times;   %In case of restart
end
