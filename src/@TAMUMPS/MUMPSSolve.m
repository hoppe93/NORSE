function f = MUMPSSolve(o, A, b)
    % Factorizes and solves the given linear system
    %    
    %   A * x = b
    %
    % where x is the unknown, using MUMPS.
    %
    % Usage:
    %    MUMPSSolve(A, b)
    %
    % A  -- (Sparse) matrix representing LHS
    % b  -- Vector representing the RHS
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    o.mumpsInfo.JOB = 5;
    o.mumpsInfo.RHS = full(b);
    o.mumpsInfo     = dmumps(o.mumpsInfo, A);
    
    f = o.mumpsInfo.SOL;

    % Check for MUMPS errors/warnings
    if o.mumpsInfo.INFOG(1) ~= 0
        if o.mumpsInfo.INFOG(1) > 0 % WARNING
            if (o.mumpsInfo.INFOG(1) & 2) == 1
                warning('MUMPS warning 2: During error analysis the max-norm of the computed solution is close to zero.');
            end

            if (o.mumpsInfo.INFOG(1) & 8) == 1
                warning(['MUMPS warning 8: Returned from iterative refinement without convergence. More than ',num2str(o.mumpsInfo.ICNTL(10)),' iterations required.']);
            end
        else    % ERROR
            switch o.mumpsInfo.INFOG(1)
                case -6
                    error(['MUMPS error -5: Matrix is singular in structure. Rank = ',num2str(o.mumpsInfo.INFOG(2))]);
                case -8
                    error('MUMPS error -8: Main internal integer workarray is too small for factorization. Increase ICNTL(14) and try again.');
                case -9
                    error('MUMPS error -9: Main internal real workarray is too small. Consult the MUMPS manual for further information.');
                case -10
                    error('MUMPS error -10: Numerically singular matrix.');
                case -11
                    error('MUMPS error -11: Internal real workarray S or LWK_USER too small for solution.');
                case -12
                    error('MUMPS error -12: Internal real workarray is too small for iterative refinement.');
                case -13
                    error('MUMPS error -13: Error in a Fortran ALLOCATE statement in MUMPS factorization. Try reducing the grid size.');
                otherwise
                    error(['MUMPS returned with an unrecognized error code: ',num2str(o.mumpsInfo.INFOG(1)),'. Please consult the MUMPS Manual for further information.']);
            end
        end
    end

end
