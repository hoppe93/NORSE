function val = Evaluate(o,t)
    % Returns the function value at the times specified by the
    % vector t.
    %
    % Usage:
    %   val = Evaluate(t)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if o.isScalar
         val = o.parameterSpline*ones(size(t));            
    else
        val = ppval(o.parameterSpline,t); %Evaluate the spline
        if o.isConstantOutsideRange
            val(t<=o.tMin) = ppval(o.parameterSpline,o.tMin);                
            val(t>=o.tMax) = ppval(o.parameterSpline,o.tMax);
        end
    end
end
