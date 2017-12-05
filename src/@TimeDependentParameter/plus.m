function val = plus(o,scal)
    % Defines the following addition operations (s is a scalar and
    % par a TimeDependentParameter):
    %   par+s=s+par & par1+par2
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isnumeric(o) && isscalar(o)
        % This corresponds to s+par. Call the function again, but
        % switch the arguments.
        val = plus(scal,o); 
    elseif isa(scal,'TimeDependentParameter')
        %Defines addition of TimeDependentParameter objects
        if scal.isScalar
            val = o+scal.parameterSpline; %Properties from o (not fact!) are inherited.                    
        elseif o.isScalar
            val = scal+o.parameterSpline; %Properties from fact (not o!) are inherited.                    
        else
            %Both are time-dependent -- create a new object defined
            %at all the time points of the two input objects
            oT = o.parameterSpline.breaks;
            scalT = scal.parameterSpline.breaks;
            times = unique(sort([oT,scalT]));
            vals = o.Evaluate(times)+scal.Evaluate(times);
            val = TimeDependentParameter(times,vals,...
                        o.isStepWiseConst,o.isConstantOutsideRange);
        end
    elseif isnumeric(scal) && isscalar(scal)    
        val = o.ApplyBinaryOperator(@plus,scal);                
    else
        error('Only addition of a scalar or another TimeDependentParameter object is defined.');
    end
end
