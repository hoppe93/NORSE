function ResetTimings(o)
    % Resets the various timers used for information about the contribution
    % to the runtime from different parts of the calculation.
    %
    % Usage:
    %   ResetTimings()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    names = fieldnames(o.timing);
    for i = 1:numel(names) %Loop through the timers and reset them
        o.timing.(names{i}) = 0;
    end    
    o.Print('Successfully reset the runtime timers.\n');
end
