function CleanUp(o)
    % Clear large internal variables that are not needed when
    % saving the NORSE object.
    % Usage:
    %   CleanUp()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    CleanUp@TimeAdvance(o);
    o.L = [];
    o.U = [];
    o.P = [];
    o.Q = [];
