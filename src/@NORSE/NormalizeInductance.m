function LHat = NormalizeInductance(o,L,a,R)
    % Normalizes the inductance L in Henry to the internal LHat used in
    % NORSE. The other physical parameters must be specified before calling
    % this function, since some of the normalization depends on them.
    %
    % Usage:
    %   LHat = NormalizeInductance(L,a,R)
    %
    % L may be a scalar, a vector, or a TimeDependentParameter. a and R are
    % the minor and major radii, respectively.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    e2 = (o.constants.e)^2;
    m = o.constants.m;
    
    if isempty(o.fM0)
        o.ProcessPhysicalParameters; 
    end
    
    LHat = L*o.fM0 * 0.5*e2*a*a/(m*R);
end
