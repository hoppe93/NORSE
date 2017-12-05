function SynchSettings(o)
    % Makes sure settings are consistent with those in NORSE.
    % Should be called before the calculation starts (or at restart).
    %
    % Usage:
    %   SynchSettings()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    oN = o.norse;
    
    %Synch the general time-advance settings
    SynchSettings@TimeAdvance(o);
        
    %Synch the settings specific to the adaptive time-advance
    o.GMRESTolerance              = oN.GMRESTolerance;
    o.nStepsBetweenFactorizations = oN.nStepsBetweenFactorizations;    
end
