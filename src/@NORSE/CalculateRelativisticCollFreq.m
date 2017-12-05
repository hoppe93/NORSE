function nu = CalculateRelativisticCollFreq(o,varargin)
    % Calculates the relativistic collision frequency nu.
    %
    % Usage:
    %   nu = CalculateRelativisticCollFreq() 
    %   nu = CalculateRelativisticCollFreq(T,n)
    %
    % The first call uses the T & n set in o, while the second call
    % uses the provided T & n. T & n may be scalars, vectors of the
    % same size or TimeDependentParameters.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Take advantage of the fact that nu is just E_c*e/mc
    e = o.constants.e;            
    m = o.constants.m;
    c = o.constants.c;                        

    switch nargin
        case 1 %Use the parameters set in the object 
            nu = o.CalculateEc() * e/(m*c);
        case 3 %Use the provided T and n
            nu = o.CalculateEc(varargin{1},varargin{2}) * e/(m*c);                    
        otherwise
            error('Invalid number of input arguments');
    end
end
