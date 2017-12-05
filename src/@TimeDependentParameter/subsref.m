function varargout = subsref(o,s)
    %Overload the standard Matlab array indexing to evaluate the
    %spline at the speficied time. This makes calls of the type 
    %par(t1) and par([t1,t2,...]) possible.
    %
    % Allowed calls:
    %   par(t) -- Evaluates the spline at the times in vector t
    %   par.MethodOrProperty()     -- Call class method or property
    %   par.MethodOrProperty(args) -- Call with arguments
    %   par.property.field     -- Access fields of a struct
    %
    % All other operations are undefined and will produce an error.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    errT = 'Invalid indexing operation. TimeDependentParameter only supports limited indexing functionality.';

    switch s(1).type
        case '()' %Evaluate the spline at the provided times
            if numel(s) == 1 && numel(s.subs) == 1 && ~isempty(s.subs{1})
                t = s.subs{1};                
                varargout{:} = o.Evaluate(t);
            else
                error(errT);
            end
        case '.' %Access class methods and properties
            switch numel(s)
                case 1 %Method/property reference without arguments
                    varargout{:} = o.(s(1).subs);
                case 2 
                    if strcmp(s(2).type,'()') %Reference with arguments to the method/property
                        varargout{:} = o.(s(1).subs)(s(2).subs{:});
                    elseif strcmp(s(2).type,'.') %Reference fields of a property that is a struct
                        varargout{:} = o.(s(1).subs).(s(2).subs);
                    else
                        error(errT); 
                    end
                otherwise
                    error(errT);
            end
        otherwise
            error(errT);
    end

    %Only provide output if the user requests it
    if numel(varargout) == 1 && isempty(varargout{1})
        varargout = {}; 
    end
end
