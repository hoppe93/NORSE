function PrintIterationWarning(o,iter)
    % Prints a warning if the GMRES made zero iterations.
    %
    % Usage:
    %   PrintConvergenceWarning(iter)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    if iter(1) == 0 && iter(2) == 0 && ~o.iterationsFlag
        warning('GMRES returned initial guess as solution.');
        o.iterationsFlag = 1;
    end
end