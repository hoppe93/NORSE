classdef TimeDependentParameter
    % TIMEDEPENDENTPARAMETER -- Class for defining time-dependent 
    %                           parameters in NORSE. 
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Description and usage: 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Allows to specify the parameter values (f) at arbitrary time points
    % (t) and storing it in a spline, which can be evaluated in NORSE as
    % neeeded. f and t must be of the same size, and may be scalar.
    %
    % Usage:
    %   par = TimeDependentParameter(t,f); 
    %   par = TimeDependentParameter(t,f,isStepWiseConst,isConstOutsideRange); 
    %
    % The class supports simple indexing and many common operations: 
    %   par = TimeDependentParameter(t,f); 
    %   val = par(t1);    %Evaluates the spline at t=t1 
    %   par2 = par*par;   %par2 is a new TimeDependentParameter such that 
    %                     %par2(tt) = par(tt)*par(tt) at any time tt. 
    %    
    % If isStepWiseConst=1 is passed, f is treated as constant in between
    % the supplied data points. The result is a step-like behavior of f 
    % (not default). This also automatically sets isConstOutsideRange=1. f
    % will change value exactly at the specified times. 
    %   Example: 
    %       t=[0,1,2], f=[0,10,0] 
    %     will result in par(0.99)=0, par(1)=10, par(1.99)=10, par(>=2)=0.
    %       
    %
    % If isConstOutsideRange=1 is passed (default), the function is assumed 
    % constant outside the specified time interval:
    %   par(t0) == par(t(1)) 
    % for t0<t(1), and 
    %   par(t2) == par(t(end))
    % for t2>t(end). Otherwise, the spline is used for extrapolation.
    %
    %
    % The following calls are allowed:
    %   par(t) -- Evaluates the spline at the times in vector t
    %   par.MethodOrProperty()     -- Call class method or property
    %   par.MethodOrProperty(args) -- Call with arguments
    %   par.property.field         -- Access fields of a property struct
    % The following operations are supported (see the respective class
    % method for details):
    %   +, -, *, /, log, log10, power (^), sqrt, double & besselk.
    %
    % All other calls and operations are undefined and will produce an
    % error, however additional unary or binary operations can easily be
    % included by calling the class methods ApplyUnaryOperator(f) and
    % ApplyBinaryOperator(f,b).
    %
    % The boolean property isScalar specifies whether the
    % timeDependentParameter describes a constant or is time dependent.
    %    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %%% Settings and interface methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties 
        isStepWiseConst = 0
        isConstantOutsideRange = 1
    end
    
    methods 
        function o = TimeDependentParameter(t,f,varargin)
            % Constructor. 
            %
            % Usage:
            %   par = TimeDependentParameter(t,f); 
            %   par = TimeDependentParameter(t,f,isStepWiseConst,isConstOutsideRange); 
            
            if ~all(size(t) == size(f)) || ~isvector(t) || ~isnumeric(t) || ~isnumeric(f)
                error('Invalid input argument dimenstions.');
            end
            if nargin >= 4  
                if isscalar(varargin{1}) && isnumeric(varargin{1}) ...
                        && isscalar(varargin{2}) && isnumeric(varargin{2}) 
                    o.isStepWiseConst = varargin{1};
                    if o.isStepWiseConst
                        o.isConstantOutsideRange = 1; 
                    else
                        o.isConstantOutsideRange = varargin{2}; 
                    end
                else
                    error('Invalid third or fourth input argument.'); 
                end
            end
            switch numel(t)
                case 0
                    o.isScalar = 1;
                    o.parameterSpline = 0;                    
                case 1
                    o.isScalar = 1;
                    o.parameterSpline = f; %Store the value directly
                    o.tMin = t;
                    o.tMax = t;
                otherwise                    
                    if o.isStepWiseConst
                        %Determine where f changes value and make a spline
                        %consisting of a series of constant polynomials 
                        %(steps).
                        breaks = find(diff(f)) + 1;
                        breaksWEdges = unique([1;breaks(:);numel(f)]);
                        fs = f(breaksWEdges);
                        o.parameterSpline = ppmak([t(breaksWEdges),t(end)+eps],fs);
                    else
                        %A continuous function -- use a spline or a pchip.
                        %The latter is better behaved since it is constant
                        %on intervals whose endpoints are the same.
                        o.parameterSpline = pchip(t,f); 
%                         o.parameterSpline = spline(t,f);                     
                    end 
                    o.tMin = t(1);
                    o.tMax = t(end);
            end
        end
        
        val = Evaluate(o,t)
            % Returns the function value at the times specified by the
            % vector t.
        o = Rescale(o,tFac)
            % Rescales the time axis so that tNew = tFac*tOld.
        out = Visualize(o,varargin)
            % Plots the spline and marks the specified data points.
        bools = HasChanged(o,ts)
            % Takes an array of times ts and determines if the parameter
            % has changed inbetween its entries.
    end
    
    %%% Internal properties and methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties 
        parameterSpline
        tMin = 0
        tMax = 0
        isScalar = 0         
    end
    
    methods %%% Utilities for indexing and overloaded operators %%% 
        varargout = subsref(o,s)
            %Overload the standard Matlab array indexing to evaluate the
            %spline at the speficied time. This makes calls of the type 
            %par(t1) and par([t1,t2,...]) possible.
        val = ApplyUnaryOperator(o,f)
            % Applies a unary operator to a TimeDependentParameter. 
        val = ApplyBinaryOperator(o,f,b)
            % Applies a binary operator to a TimeDependentParameter. 
    end
    
    methods %%% Overloaded operators %%%
        % The operations +, -, * and / are defined below.
        %
        % When operating on two TimeDependentParameters, the operation is
        % performed on the function values at all time points (for any
        % object) where the spline is specified. The result is a new
        % TimeDependentParameter, defined at all time points of either
        % initial object.
        %
        % Note that, in general, the operations do not commute since the
        % properties isStepWiseConst and isConstantOutsideRange are
        % inherited from the first non-constant object.
        val = mrdivide(o,denom)
            % Defines the following division operations (s is a scalar and
            % par a TimeDependentParameter):
            %   par/s, s/par & par/par
        val = mtimes(o,fact)
            % Defines the following multiplication operations:
            %   par*s=s*par & par1*par2
        val = plus(o,scal)
            % Defines the following addition operations:
            %   par+s=s+par & par1+par2
        val = minus(o,scal)
            % Defines the following subtraction operations:
            %   par-s, s-par & par1-par2
        
        
        %Other unary or binary operators -- perform the corresponding
        %operation on the points where the spline is defined.
        function val = log(o)
            val = o.ApplyUnaryOperator(@log);
        end
        
        function val = log10(o)
            val = o.ApplyUnaryOperator(@log10);
        end        
        
        function val = power(o,xp)
            %Only supports scalar exponents
            
            if ~isa(o,'TimeDependentParameter') || ~isnumeric(xp) || ~isscalar(xp)
                error('The operation is not defined.');
            end
            val = o.ApplyBinaryOperator(@power,xp);
        end
        
        function val = sqrt(o)
            val = o.ApplyUnaryOperator(@sqrt);
        end        
        
        function val = double(o)
            %Only works on scalar objects
            
            if o.isScalar
                val = o.parameterSpline;
            else
                error('Conversion to double is only possible for scalar TimeDependentParameter objects.');
            end
        end
        
        function val = besselk(nu,o,varargin)
            if nargin == 3
                f = @(x) besselk(nu,x,varargin{1});
            else
                f = @(x) besselk(nu,x);
            end
            val = o.ApplyUnaryOperator(f);
        end
    end    
end
