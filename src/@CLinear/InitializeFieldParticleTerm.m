function InitializeFieldParticleTerm(o) 
    % Calculates the time-independent prefactors for the field-particle
    % term of the electron-electron collision operator.
    %
    % Usage:
    %   InitializeFieldParticleTerm()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%% Preliminaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    oN = o.norse;
    oG = o.grid;
    p     = oG.pBig;
    g = oG.gammaBig;
    xi    = oG.xiBig;

    %Helper function for diagonal matrices, to reduce clutter
    dia = @(v) spdiags(v,0,oG.matSize,oG.matSize); 

    %Miscellaneous grid quantities
    p2     = p.*p;
    p3     = p.*p2;
    g2     = 1+p2;
    divP   = 1./p; 
    divG   = 1./g;
    divP2  = divP.*divP;                                    
    p2DivG = p2.*divG;
    divGP  = divP.*divG;  
    oMXiSq = 1-xi.*xi;
    oMXiSqDivP2 = oMXiSq.*divP2;


    %%% Calculate the prefactors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This will multiply the corresponding potential weighted with 
    % exp((1-gamma)/Theta). This is because it is numerically preferrable
    % to apply the exponential factor in the Legendre basis since it avoids
    % the creation of large full matrices than then need to be multiplied
    % by the exponential factor (and most of the entries thrown away).
    
    %Calculate the normalization factor 
    o.preFactorFP = oN.nBar*oN.lnLambdaBar/(oN.Theta*oN.kappa); 

    %%% Upsilon0
    ddp  = dia(p).*oG.ddpMat;     
    noD1 = 3*speye(oG.matSize);
    noD2 = dia(p2DivG);
    
    o.dUpsilon0.invTh  = ddp + noD1; %This should be multiplied by 1/Theta(t)
    o.dUpsilon0.invTh2 = noD2; %This should be multiplied by 1/Theta(t)^2
    
    %%% Upsilon2
    ddp  = -dia(8*p)*oG.ddpMat;    
    noD1 = -24*speye(oG.matSize);
    noD2 = -dia(8*p2DivG);
    
    o.dUpsilon2.invTh  = ddp + noD1;
    o.dUpsilon2.invTh2 = noD2;
    
    %%% UpsilonMin
    d2dp2_1 = dia(2*g2)*oG.d2dp2Mat;
    d2dp2_2 = dia(g.*p2)*oG.d2dp2Mat;     
    d2dxi2  = dia(2*oMXiSqDivP2)*oG.d2dxi2Mat; 
    ddp1    = dia(4*divP+8*p)*oG.ddpMat;
    ddp2    = dia(p.*divG.*(4+5*p2))*oG.ddpMat;
    ddp3    = dia(2*p3)*oG.ddpMat;    
    ddxi    = -dia(4*xi.*divP2)*oG.ddxiMat; 
    noD2    = dia(divG.*(6+8*p2));
    noD3    = dia(p2DivG.*divG.*(3+2*p2+g.*p));
    noD4    = dia(p2DivG.*p2DivG);
    
    o.dUpsilonMin.invTh  = d2dp2_1 + d2dxi2 + ddp1 + ddxi;
    o.dUpsilonMin.invTh2 = d2dp2_2 + ddp2 + noD2;
    o.dUpsilonMin.invTh3 = ddp3 + noD3;
    o.dUpsilonMin.invTh4 = noD4;
    
    %%% UpsilonPlus
    o.dUpsilonPlus = -dia(p2DivG); %This should be multiplied by 1/Theta(t)^2
    
    %%% Pi
    d2dp2  = -dia(g)*oG.d2dp2Mat; 
    d2dxi2 = -dia(oMXiSqDivP2.*divG)*oG.d2dxi2Mat; 
    ddp0   = -dia(divGP.*(2+3*p2))*oG.ddpMat; 
    ddp1   = -dia(p)*oG.ddpMat;
    ddxi   = dia(2*xi.*divGP.*divP)*oG.ddxiMat; 
    noD1   = -dia(divG.*divG.*(3+3*p2));    
    
    o.dPi.noTh   = d2dp2 + d2dxi2 + ddp0 + ddxi; %This has no Theta dependence
    o.dPi.invTh  = ddp1 + noD1;
    
    %%% Calculate matrices for the combined potentials %%%%%%%%%%%%%%%%%%%%
    oP = o.potentials;
    o.UpsilonMinMatrix  = 4*oP.Upsilon2Matrix - oP.Upsilon1Matrix;                   
    o.UpsilonPlusMatrix = 4*oP.Upsilon2Matrix + oP.Upsilon1Matrix;
end
