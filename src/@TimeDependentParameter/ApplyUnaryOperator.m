function val = ApplyUnaryOperator(o,f)
    % Applies a unary operator to a TimeDependentParameter. f must
    % be a function handle to a unary operator such as log or sin.
    %
    % Usage:
    %   val = ApplyUnaryOperator(f)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if o.isScalar
        val = o;
        val.parameterSpline = f(val.parameterSpline);
    elseif o.isStepWiseConst %Change value at defined time points
        val = o;
        val.parameterSpline.coefs = f(val.parameterSpline.coefs);
    else %Smooth function -- make a new object
        breaks = o.parameterSpline.breaks;
        val = TimeDependentParameter(breaks,f(o.Evaluate(breaks)),...
                        o.isStepWiseConst,o.isConstantOutsideRange);
    end
end
