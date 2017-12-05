function f = PreconditionerFunc(o,rhs)
    % Efficiently computes the preconditioner used in the iterative
    % time advance schemes. Significantly speeds up the computation,
    % compared to supplying the L & U matrices directly.
    %
    % Usage: 
    %   f = PreconditionerFunc(rhs)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    f = o.Q*(o.U\(o.L\(o.P*rhs)));
end    
