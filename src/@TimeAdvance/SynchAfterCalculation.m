function SynchAfterCalculation(o)
    % Updates some fields in the NORSE object to reflect the final state of
    % the time advance.
    %
    % Usage: 
    %   SynchAfterCalculation()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    oN            = o.norse;
    
    oN.times      = o.times;
    oN.nTimeSteps = o.nTimeSteps;
    oN.nSaveSteps = o.nSaveSteps;    
end
