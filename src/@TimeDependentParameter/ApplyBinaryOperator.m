function val = ApplyBinaryOperator(o,f,b)
    % Applies a binary operator to a TimeDependentParameter. f must
    % be a function handle to a binary operator (such as power)
    % which takes a second argument b.
    %
    % Usage:
    %   val = ApplyBinaryOperator(f,b)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if o.isScalar
        val = o;
        val.parameterSpline = f(val.parameterSpline,b);
    elseif o.isStepWiseConst %Change value at defined time points
        val = o;
        val.parameterSpline.coefs = f(val.parameterSpline.coefs,b);
    else %Smooth function -- make a new object
        breaks = o.parameterSpline.breaks;
        val = TimeDependentParameter(breaks,f(o.Evaluate(breaks),b),...
                        o.isStepWiseConst,o.isConstantOutsideRange);
    end
end        
