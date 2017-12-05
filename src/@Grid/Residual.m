function r = Residual(v,d,nXi,id)
    %Returns the residual for the optimization to find a sigmoid
    %grid mapping with a point at xi=0.
    %
    % Usage:
    %   r = Residual(v,d,nXi,id)
    %
    % r is the distance of the closest grid point id (established
    % in FindOptimalSigmoidMapping) to 0. v = [c,h], with c,d and
    % h parameters in the sigmoid mapping.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    c = v(1); 
    h = v(2);

    %Generate a uniform grid with the given c and h values
    lowerSBound = log((h-1)/(c-h+1))/d;
    upperSBound = log((h+1)/(c-h-1))/d;            
    s = linspace(lowerSBound,upperSBound,nXi);

    %Do the sigmoid mapping. The residual is the distance to 0 of
    %the closest point.
    r = c./(1+exp(-d*s(id)))-h;
    if isreal(r)
        r = abs(r);
    else
        r=Inf;
    end
end
