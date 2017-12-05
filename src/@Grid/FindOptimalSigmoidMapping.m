function [c,h] = FindOptimalSigmoidMapping(c,d,h,nXi)
    %Finds a grid mapping with a point at xi=0 and endpoints at
    %xi=-1,1, given an initial guess for the parameters c,d and h
    %in the sigmoid mapping for the xi grid. Returns the
    %corresponding values for c and h.
    %
    % Usage:
    %   [c,h] = FindOptimalGridMapping(c,d,h,nXi)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Create a uniform grid with endpoints appropriate for the
    %sigmoid mapping
    lowerSBound = log((h-1)/(c-h+1))/d;
    upperSBound = log((h+1)/(c-h-1))/d;            
    s = linspace(lowerSBound,upperSBound,nXi);

    %Find the grid point closest to 0. This is the point in the
    %grid we will use for the optimization of c and h.
    [~,id] = min(abs(c./(1+exp(-d*s))-h)); 

    fToMinimize = @(v) Grid.Residual(v,d,nXi,id);
    op = optimset('TolFun',1e-14); 

    v = fminsearch(fToMinimize,[c,h],op);

    c = v(1);
    h = v(2);
end
