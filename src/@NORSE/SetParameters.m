function SetParameters(o,nP,nXi,nL,pMax,dt,tMax,T,n,Z,EHat,varargin)
    %This function is handy for setting several parameters in one
    %call. The physical parameters can be scalars or
    %TimeDependentParameter objects.
    %
    %Usage:
    %   SetParameters(nP,nXi,nL,pMax,dt,tMax,T,n,Z,EHat)
    %   SetParameters(nP,nXi,nL,pMax,dt,tMax,T,n,Z,EHat,B)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    o.nP   = nP;
    o.nXi  = nXi;
    o.nL   = nL;
    o.pMax = pMax;
    o.dt   = dt;
    o.tMax = tMax;            
    o.T    = T;
    o.n    = n;
    o.Z    = Z;
    o.EHat = EHat;
    if nargin >= 12
        o.B = varargin{1};
    end
end
