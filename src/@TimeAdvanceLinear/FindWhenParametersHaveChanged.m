function parsHaveChanged = FindWhenParametersHaveChanged(o,firstStep)
    % Evaluates when the time-dependent physical parameters have changed so
    % that the matrix will need to be rebuilt. The result is returned as a
    % struct with fields for each parameter consisting of a logical array
    % with an entry for each time step.
    %
    % Usage:
    %   parsHaveChanged = FindWhenParametersHaveChange(firstStep)    
    %    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    oN = o.norse;    
    ts = o.allTimes(firstStep-1:end);    
    parsHaveChanged     = struct();
    
    parsHaveChanged.T   = oN.Theta.HasChanged(ts);
    parsHaveChanged.n   = oN.nBar.HasChanged(ts);
    parsHaveChanged.E   = oN.EHat.HasChanged(ts);
    parsHaveChanged.Z   = oN.Z.HasChanged(ts);
    
    parsHaveChanged.any = parsHaveChanged.T | parsHaveChanged.n | ...
                          parsHaveChanged.E | parsHaveChanged.Z;
    parsHaveChanged.any(firstStep) = true;
                                %Make sure the matrix will be built in the
                                %first time step
end
