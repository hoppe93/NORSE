function PrintConvergenceWarning(o,flag)
    % Prints a warning if the GMRES iteration raises an error flag.
    %
    % Usage:
    %   PrintConvergenceWarning(gmresFlag)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if flag && ~o.convergenceFlag
        switch flag
            case 1
                warning('GMRES failed to converge (error flag: %d)',flag) 
            case 2
                warning('The preconditioner to GMRES was illconditioned (error flag: %d)',flag) 
            case 3
                warning('GMRES stagnated (error flag: %d)',flag) 
        end
        o.convergenceFlag = 1;
    end
end
