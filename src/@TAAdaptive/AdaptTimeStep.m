function dt = AdaptTimeStep(o,dt,nIterations)
    % Adapts time step based on the number of GMRES steps to
    % convergence. 
    %
    % Usage:
    %   dt = AdaptTimeStep(dt,nGMRESIterations)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Set the desired number of iterations
    optim = 6; 

    if o.dtIncreaseLimit
        %Restrict the maximum allowed time step to a specified
        %multiple of the initial time step
        isTooLarge = dt > o.dtIncreaseLimit*o.dt;
    else
        isTooLarge = false;
    end

    %Adapt the time step
    if ~nIterations
        %First time step, nothing to do
    elseif nIterations > optim+9
        dt = 0.4*dt;
    elseif nIterations > optim+4
        dt = 0.65*dt;
    elseif nIterations > optim+1
        dt = 0.9*dt;
    elseif nIterations < optim-3 && ~isTooLarge
        dt = 2*dt;
    elseif nIterations < optim-2 && ~isTooLarge
        dt = 1.5*dt;
    elseif nIterations < optim-1 && ~isTooLarge
        dt = 1.25*dt;                        
    end
end
