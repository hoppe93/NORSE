function val = mrdivide(o,denom)
    % Defines the following division operations (s is a scalar and
    % par a TimeDependentParameter):
    %   par/s, s/par & par/par
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isnumeric(o) && isscalar(o) %s/par
        %We have s/par, since we know that at least one of the
        %arguments must be the object.
        val = denom;
        if denom.isScalar                    
            val.parameterSpline = o./val.parameterSpline;  
        elseif denom.isStepWiseConst                    
            val.parameterSpline.coefs = o./val.parameterSpline.coefs;  
        else %Continuous function --  create a new object
            breaks = val.parameterSpline.breaks;
            val = TimeDependentParameter(breaks,o./val.Evaluate(breaks),...
                        val.isStepWiseConst,val.isConstantOutsideRange);                    
        end
    elseif isa(denom,'TimeDependentParameter') %par/par
        %Defines division between TimeDependentParameter objects
        if denom.isScalar
            val = o/denom.parameterSpline; %Properties from o (not denom!) are inherited.                    
        elseif o.isScalar
            val = o.parameterSpline/denom; %Properties from denom (not o!) are inherited.                    
        else
            %Both are time-dependent -- create a new object defined
            %at all the time points of the two input objects
            oT = o.parameterSpline.breaks;
            denomT = denom.parameterSpline.breaks;
            times = unique(sort([oT,denomT]));
            vals = o.Evaluate(times)./denom.Evaluate(times);
            val = TimeDependentParameter(times,vals,...
                        o.isStepWiseConst,o.isConstantOutsideRange);
        end
    elseif isnumeric(denom) && isscalar(denom) %par/s
        val = o.ApplyBinaryOperator(@mrdivide,denom);                
    else
        error('Only division involving a scalar or another TimeDependentParameter object is defined.');                
    end
end
