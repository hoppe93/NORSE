function Update(o,f,t,tOld)
    % Updates and returns the properly scaled sink.
    %
    % Usage:
    %   Update(f,t,tOld)
    %
    % f and t are the distribution and time at the current time
    % step. tOld is the time at the previous time step.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    o.CalculateEnergyChangeMagnitude(f,t,tOld);
    o.sink = o.energyChangeMagnitude.sink*o.heatSinkOperator;
end
