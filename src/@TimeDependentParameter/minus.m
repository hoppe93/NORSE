function val = minus(o,scal)
    % Defines the following subtraction operations (s is a scalar
    % and par a TimeDependentParameter):
    %   par-s, s-par & par1-par2
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isa(o,'TimeDependentParameter') && isa(scal,'TimeDependentParameter')
        %par-par; defines subtraction of TimeDependentParameter objects
        if scal.isScalar
            val = o-scal.parameterSpline; %Properties from o (not fact!) are inherited.                    
        elseif o.isScalar
            val = o.parameterSpline-scal; %Properties from fact (not o!) are inherited.                    
        else
            %Both are time-dependent -- create a new object defined
            %at all the time points of the two input objects
            oT = o.parameterSpline.breaks;
            scalT = scal.parameterSpline.breaks;
            times = unique(sort([oT,scalT]));
            vals = o.Evaluate(times)-scal.Evaluate(times);
            val = TimeDependentParameter(times,vals,...
                        o.isStepWiseConst,o.isConstantOutsideRange);
        end
    elseif isa(o,'TimeDependentParameter')
        %par-s
        if isnumeric(scal) && isscalar(scal)
            val = o.ApplyBinaryOperator(@minus,scal);
        else
            error('This subtraction operation is not defined.');
        end
    else
        %s-par
        if isnumeric(o) && isscalar(o)                    
            val = scal;
            if val.isScalar
                val.parameterSpline = o-val.parameterSpline;
            elseif val.isStepWiseConst
                val.parameterSpline.coefs = o-val.parameterSpline.coefs;
            else
                breaks = val.parameterSpline.breaks;
                val = TimeDependentParameter(breaks,o-val.Evaluate(breaks),...
                        val.isStepWiseConst,val.isConstantOutsideRange);
            end
        else
            error('This subtraction operation is not defined.');
        end
    end
end
