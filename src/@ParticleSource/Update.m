function Update(o,t,tOld,fOld,matrix,rhs)
    % Updates and returns the properly scaled sink.
    %
    % Usage:
    %   Update(f,t,tOld)
    %
    % f and t are the distribution and time at the current time
    % step. tOld is the time at the previous time step.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    o.BuildAndCalculateMagnitude(t,tOld,fOld,matrix,rhs);
    o.source = o.densityChangeMagnitude*o.particleSourceOperator;
end
