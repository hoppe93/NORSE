function CalculateSpecialFunctions(o,t) 
    % Determines the special functions G_K, G_Para, G_Perp (as well as the
    % derivatives dG_K/dp and dG_Para/dp) needed for the test-particle term
    % of the linear e-e collision operator at time t.
    %
    % Usage:
    %   CalculateSpecialFunctions(t)
    %    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    oN = o.norse;
    oG = o.grid;
    th = oN.Theta(t);
    ka = oN.kappa(t);
    
    %%% Evaluate the special functions on the 1D p grid %%%%%%%%%%%%%%%%%%%
    p  = oG.p;
    p2 = p.*p;
    g  = oG.gamma; 
    g2 = g.*g;
    
    % Evaluate psi0 and psi1 efficiently or accurately. The quick
    % evaluation is orders of magnitude faster
    if o.specialFunctionEvaluation
        %Evaluate psi0 and psi1 using Simpson's rule -- quick
        psi0F = @(p) (exp((1-sqrt(1+p.^2))/th)./sqrt(1+p.^2));
        psi1F = @(p)  exp((1-sqrt(1+p.^2))/th);           
        psi0  = o.CumSimpsonsRule(psi0F,p);
        psi1  = o.CumSimpsonsRule(psi1F,p);
    else
        %This implementation is based on the work of Eero Hirvijoki and 
        %uses MATLAB's integral function (i.e. it will be correct to any 
        %desired accuracy) -- slow
        
        %Evaluate psi0 and psi1 using quad
        psi0Int = @(x) exp((1-sqrt(1+x.^2))/th)./sqrt(1+x.^2);
        psi1F   = @(x) exp((1-sqrt(1+x.^2))/th);
        
        % calculate the special functions
        psi0    = zeros(size(p));
        psi1    = zeros(size(p));    
        as      = [0;p(1:end-1)];
        bs      = p;
        psi0(1) = integral( @(x) psi0Int(x),as(1),bs(1) );
        psi1(1) = integral( @(x) psi1F(x),as(1),bs(1) );
        
        for i=2:numel(p)        
            psi0(i) = integral( @(x) psi0Int(x),as(i),bs(i) ) + psi0(i-1);
            psi1(i) = integral( @(x) psi1F(x),  as(i),bs(i) ) + psi1(i-1);
        end  
    end
        
    expFactor = exp((1-g)/th);
    GK    = ( g2.*psi1 - th*psi0 + (th*g-1).*p.*expFactor )./(p2*ka);
    GPara = th*g.*GK./p;
    GPerp = ( (p2.*g2+th^2).*psi0 + th*(2*p2.*p2-1).*psi1 ...
                + th*p.*g.*(1+th*(2*p2-1)).*expFactor          )...
              ./(2*g.*p2.*p*ka);
          
    dGKDP = 1/ka *(  2./(p.*p2).*(th*psi0-psi1) ...
                   + 1./p2.*(2-2*th./g + p2./(g*th)).*expFactor);
    dGParaDP = -th./(g.*p2) .* GK  +  th*g./p .* dGKDP;
    
%     % set the values for p=0 separately  
%     GK(1)    = 0;
%     GPara(1) = (1+2*th+2*th^2)/(3*ka);
%     GPerp(1) = (4+5*th+8*th^2)/(12*ka);
%     dGKDP(1) = 0; %Do we need to calculate these?
%     dGParaDP(1) = 0; %Do we need calculate these?
    
    %%% Map the results to a vector on the 2D grid %%%%%%%%%%%%%%%%%%%%%%%%
    onesVec    = ones(1,oN.nXi);
    o.GK       = oG.BuildBigVector(onesVec,GK);
    o.GPara    = oG.BuildBigVector(onesVec,GPara);
    o.GPerp    = oG.BuildBigVector(onesVec,GPerp);
    o.dGKDP    = oG.BuildBigVector(onesVec,dGKDP);
    o.dGParaDP = oG.BuildBigVector(onesVec,dGParaDP);
end
