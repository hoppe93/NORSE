function MUMPSAnalyze(o, matrix)
    % Run the MUMPS analysis stage. This function checks for errors in
    % MUMPS and outputs relevant warnings/error messages.
    %
    % Usage:
    %    MUMPSAnalyze(matrix)
    %
    % matrix -- (Sparse) matrix to analyze
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    o.mumpsInfo.JOB = 1;
    o.mumpsInfo     = dmumps(o.mumpsInfo, matrix);

    if o.mumpsInfo.INFOG(1) ~= 0
        switch o.mumpsInfo.INFOG(1)
            case -5
                error(['MUMPS error -5: Problem of real workspace allocation of size ',num2str(o.mumpsInfo.INFOG(2)),' during analysis.']);
            otherwise
                error(['MUMPS returned with an unrecognized error code: ',num2str(o.mumpsInfo.INFOG(1)),'.']);
        end
    end

end
