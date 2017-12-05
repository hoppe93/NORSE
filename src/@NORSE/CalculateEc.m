function Ec = CalculateEc(o,varargin)
    % Calculates the critical field of Connor & Hastie [NF, 15
    % (1975) 415] in V/m.
    %
    % Usage:
    %   Ec = CalculateEc()    -- uses the T & n set in o
    %   Ec = CalculateEc(T,n) -- uses the provided T & n
    %
    % T and n may be scalars, vectors of the same size or
    % TimeDependentParameters.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    e3 = o.constants.e^3;
    eps02 = o.constants.eps0^2;
    m = o.constants.m;
    c2 = (o.constants.c)^2;

    switch nargin
        case 1 %Use the parameters set in the object
            lnLambda = o.CalculateLnLambda();
            Ec =  o.n*lnLambda*e3 / (4*pi*eps02*m*c2);
        case 3 %Use the provided T and n
            T = varargin{1};
            n = varargin{2};
            if size(T) ~= size(n)
               error('The dimensions of the inputs T and n are inconsistent');
            end
            isArrayValued = (numel(T)>1);                         

            lnLambda = o.CalculateLnLambda(T,n);
            if isArrayValued
                %For arrays we need to use .*
                Ec =  n.*lnLambda*e3 / (4*pi*eps02*m*c2);
            else                        
                %For TimeDependentParameters we cannot use .* (also
                %works on scalars, of course)
                Ec =  n*lnLambda*e3 / (4*pi*eps02*m*c2);
            end
        otherwise
            error('Invalid number of input arguments');
    end
end 
