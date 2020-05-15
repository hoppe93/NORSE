function InitializeMUMPS(o,mumpsPath)
    % Initializes the MUMPS solver.
    %
    % Usage:
    %    InitializeMUMPS(mumpsPath)  -- Specifies the path to MUMPS.
    %                                   If empty, the path must already
    %                                   be in the MATLAB path.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % MUMPS path given?
    if ~isempty(mumpsPath)
        addpath(mumpsPath);
    else
        % Do not use MUMPS
        return
    end

    % Try to find MUMPS routines
    if exist('dmumpsmex', 'file') ~= 3
        warning(['Unable to find the MUMPS MEX file (dmumpsmex) in the specified directory: ',mumpsPath]);
        return
    end

    o.mumpsInfo = initmumps;
    o.mumpsInfo = dmumps(o.mumpsInfo);

    % Disable all output except errors
    o.mumpsInfo.ICNTL(1) = 6;   % Errors (6 = standard output)
    o.mumpsInfo.ICNTL(2) = 0;   % Diagnostics/statistics/warnings (0 = disabled)
    o.mumpsInfo.ICNTL(3) = 0;   % Global info (0 = disabled)
    o.mumpsInfo.ICNTL(4) = 1;   % Level of output (1 = errors only)
    
    o.useMumps = 1;
    disp('Using MUMPS for matrix inversions');

end
