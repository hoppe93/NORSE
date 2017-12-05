function CalculateNumberOfSteps(o)
    %Perform some initialization needed to set up the rest of the NORSE
    %calculation. Most of the actual time-advance-specific initialization
    %is done in InitializeTimeStepArrays()
    %
    % Overloads the method of the TimeAdvance superclass
   
    
    oN = o.norse;    
    if isempty(oN.nTimeSteps) && ~oN.nSaveSteps
        %Have a guess at the number of steps to save if the
        %user hasn't specified anything (and it's not a restart)
        oN.nSaveSteps = min([101,round(oN.tMax/oN.dt)+1]);
    end
end
