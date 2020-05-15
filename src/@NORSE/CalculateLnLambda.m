function lnLambda = CalculateLnLambda(o,varargin) 
    %Calculates the Coulomb logarithm. The formula used is taken
    %from the NRL formulary.
    %
    %Usage:
    %   lnLambda = CalculateLnLambda()    -- uses the n & T set in o
    %   lnLambda = CalculateLnLambda(T,n) -- uses the provided n & T
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   isArrayValued = 0;
   switch nargin
       case 1                   
           T = o.T;
           n = o.n;
       case 3
           T = varargin{1};
           n = varargin{2};
           if size(T) ~= size(n)
               error('The dimensions of the inputs T and n are inconsistent');
           end
           isArrayValued = (numel(T)>1);                                         
       otherwise
           error('Wrong number of input arguments')
   end            
    if isArrayValued
        %T and n are arrays (in general), and we need to use .*
        %and ./
        %lnLambda = 23.5-log(sqrt(1e-6*n)./(T).^(5/4))...
        %     - sqrt(1e-5 + (log(T)-2).^2/16);
        lnLambda = 14.9 - 0.5*log(n/1e20) + log(T/1e3);
    else
        %T and n are timeDependentParameters, which do not
        %support .* and ./
        %lnLambda = 23.5-log(sqrt(1e-6*n)/(T).^(5/4))...
        %     - sqrt(1e-5 + (log(T)-2).^2/16);
        lnLambda = 14.9 - 0.5*log(n/1e20) + log(T/1e3);
    end            
end
