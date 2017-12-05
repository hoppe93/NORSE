function val = mtimes(o,fact)
    % Defines the following multiplication operations (s is a
    % scalar and par a TimeDependentParameter):
    %   par*s=s*par & par1*par2
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isa(o,'TimeDependentParameter')
        % This should correspond to s*par. Call the function again,
        % but switch the arguments.
        val = mtimes(fact,o);                 
    elseif isa(fact,'TimeDependentParameter')
        %Defines multiplication of TimeDependentParameter objects
        if fact.isScalar
            val = o*fact.parameterSpline; %Properties from o (not fact!) are inherited.                    
        elseif o.isScalar
            val = fact*o.parameterSpline; %Properties from fact (not o!) are inherited.                    
        else
            %Both are time-dependent -- create a new object defined
            %at all the time points of the two input objects
            oT = o.parameterSpline.breaks;
            factT = fact.parameterSpline.breaks;
            times = unique(sort([oT,factT]));
            vals = o.Evaluate(times).*fact.Evaluate(times);
            val = TimeDependentParameter(times,vals,...
                        o.isStepWiseConst,o.isConstantOutsideRange);
        end
    elseif isnumeric(fact) && isscalar(fact)    
        val = o.ApplyBinaryOperator(@mtimes,fact);                
    else
        error('Only multiplication by a scalar or another TimeDependentParameter object is defined.');                
    end
end
