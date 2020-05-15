function BuildSourceMatrix(o)
    % Builds the avalanche source matrix on the entire grid.
    %
    % Usage:
    %   BuildSourceMatrix()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    oN = o.norse;
    
    lnLambda0 = oN.CalculateLnLambda();
    o.lnLambda0 = lnLambda0(0);
    
    % Momentum of the incoming runaway
    p = o.grid.pBig;
    xi = o.grid.xiBig;
    xi2 = xi.^2;
    p1 = 2*p.*xi ./ (1+xi2 - sqrt(p.^2+1).*(1-xi2));
    p12 = p1.^2;
    
    o.intOperator = sparse(o.grid.matSize, o.grid.matSize);
    
    % Better implementation
    gmm = o.grid.gamma;
    gammaPrimaryMax = gmm(end); % Ok as long as p(k) < p(k+1) (grid increasing)
    ximin = sqrt((gammaPrimaryMax+1)/(gammaPrimaryMax-1) * (gmm-1)./(gmm+1));
    ximax = sqrt(gmm ./ (gmm+1));
    
    % Store indices of all p(1) (minus one),
    % corresponding to different xi values
    allxi = ((1:o.grid.nXi) - 1) * (o.grid.nP-1);
    
    % Calculate the Moller factor
    moller = @(g,g1) g1.^2 ./ ((g1.^2-1).*(g-1).^2.*(g1-g).^2) .* ...
        ((g1-1).^2 - (g-1).*(g1-g)./g1.^2 .* (2*g1.^2+2*g1-1-(g-1).*(g1-g)));
    
    % Create prefactor
    gammaOne = sqrt(p12 + 1);
    prefac = 1/(2*o.lnLambda0) .* p12 .* p12 ./ ...
        (o.grid.pBig.*o.grid.gammaBig.*o.grid.xiBig);
    pf = (prefac .* moller(o.grid.gammaBig, gammaOne));
    
    % We start from ip = 2, since
    % p(1) = p_0 = 0.
    for ip=2:o.grid.nP-1
        % Index of p(ip) in o.grid.pBig
        realPIndex = ip-1;
        
        iXiSecondary = find((o.grid.xi >= ximin(ip)) & (o.grid.xi <= ximax(ip)));
        if isempty(iXiSecondary)
            continue;
        end
        
        for ixi=iXiSecondary(1):iXiSecondary(end)
            sourceIndex = (ixi-1)*(o.grid.nP-1) + realPIndex;
            p1val = p1(sourceIndex);
             
            % Find first 'p(i) > p1'
            indx = find(o.grid.p > p1val, 1) - 1;
            % Make sure p1 is on the p-grid
            if isempty(indx)
                continue;
            end
            
            % Make sure p1 is within the grid
            % (and not = p_0 = 0 (which is at p(end)))
            if (indx+1 <= o.grid.nP-1)
                c = (p1val - p(indx))/(p(indx+1)-p(indx));
                pfFactor = pf(sourceIndex);
                o.intOperator(sourceIndex, allxi + indx) = c * pfFactor .* o.grid.xiWeights;
                o.intOperator(sourceIndex, allxi + indx+1) = (1-c) * pfFactor .* o.grid.xiWeights;
            else
                o.intOperator(sourceIndex, allxi + indx) = pfFactor .* o.grid.xiWeights;
            end
        end
    end

end