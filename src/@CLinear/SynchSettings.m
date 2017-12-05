function SynchSettings(o)
    % Makes sure settings are consistent with those in NORSE.
    % Should be called before the calculation starts (or at restart).
    %
    % Usage:
    %   SynchSettings()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    oN = o.norse;
    
    %Synch collision-operator settings from the NORSE object
    o.potentialCutOff           = oN.potentialCutOff;
    o.specialFunctionEvaluation = oN.specialFunctionEvaluation;
    switch oN.collisionOperatorMode
        case 1
            o.includeTPTerm = 1;
            o.includeFPTerm = 1;
        case 2
            o.includeTPTerm = 1;
            o.includeFPTerm = 0;
        otherwise
            %Something has gone wrong in the object initialization
            error('Invalid collision operator mode for CLinear.')
    end
end
