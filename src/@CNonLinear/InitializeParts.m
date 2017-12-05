function InitializeParts(o) 
    % Calculates the prefactors that multiply the different
    % potentials in the various terms (grouped by the derivative of
    % f) in the collision operator. 
    %
    % Each prefactor is saved as a field in a struct -- one struct
    % for each derivative of f (class properties d2fdp2, dfdp,
    % d2fdxi2, dfdxi, d2fdpdxi, f, ionD2fdxi2 & ionDfdxi) with one
    % field for each (non-zero) potential prefactor. For instance,
    % the prefactor multiplying the potential Upsilon_{min] in the
    % dfdp term of the electron-electron collision operator is
    % saved in dfdp.UpsMin.
    %
    % Usage:
    %   InitializeParts()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Preliminaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    oN = o.norse;
    oG = o.grid;
    p     = oG.pBig;
    gamma = oG.gammaBig;
    xi    = oG.xiBig;

    %Helper function for diagonal matrices, to reduce clutter
    dia = @(v) spdiags(v,0,oG.matSize,oG.matSize); 

    %Miscellaneous grid quantities
    divP  = 1./p; 
    divG = 1./gamma;
    divP2 = divP.*divP;                                    
    pDivG = p.*divG;
    divGP = divP.*divG;
    divGP2 = divGP.*divP;
    divGP3 = divGP2.*divP;
    divGP4 = divGP3.*divP;
    gDivP = gamma.*divP;
    gDivP2 = gDivP.*divP;
    gDivP3 = gDivP2.*divP;
    g3DivP = gamma.^3.*divP;            
    oMXiSq  = 1-xi.*xi;
    oMXiSqDivP2 = oMXiSq.*divP2;
    oMXiSqDivGP2 = oMXiSqDivP2.*divG;            
    oMXiSqDivGP4 = oMXiSqDivGP2.*divP2;


    %%% Calculate the prefactors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Prefactors for the electron-electron operator

    %Calculate the constant normalization factor
    invTh0 = 1/oN.Theta(oN.referenceTime);           
    o.preFactorEE = oN.lnLambdaBar*invTh0 / besselk(2,invTh0,1);

    %d2f/dp2 term
    d2fdp2.Ups0 = dia(-gamma);
    d2fdp2.Ups2 = dia(8*gamma);
    d2fdp2.UpsMin = -dia(2*g3DivP)*oG.ddpMat ...
                    -dia(gDivP2.*oMXiSq)*oG.d2dxi2Mat ...
                    +dia(2*gDivP2.*xi)*oG.ddxiMat; 

    %df/dp term
    dfdp.Ups0 = -dia(2*divGP+3*pDivG)-dia(gamma)*oG.ddpMat; 
    dfdp.Ups1 = dia(6*gamma)*oG.ddpMat;
    dfdp.Ups2 = 8*dia(2*divGP+3*pDivG)-dia(16*gamma)*oG.ddpMat;            
    dfdp.UpsMin = dia(-2*g3DivP)*oG.d2dp2Mat...
                  -dia(2*g3DivP.*divP)*oG.ddpMat...
                  +dia(2*xi.*(2*divGP+divGP3))*oG.ddxiMat... 
                  -dia(oMXiSq.*(2*divGP+divGP3))*oG.d2dxi2Mat;
    dfdp.Pi = -dia(gamma)*oG.ddpMat;

    %d2f/dxi2 term
    d2fdxi2.UpsPlus = dia(-oMXiSqDivGP2);
    d2fdxi2.UpsMin = dia(oMXiSq.*gDivP3)*oG.ddpMat ...
                     +dia(oMXiSq.*oMXiSqDivGP4)*oG.d2dxi2Mat ...
                     -dia(xi.*oMXiSqDivGP4)*oG.ddxiMat;            

    %df/dxi term
    dfdxi.Ups0 = dia(-oMXiSqDivGP2)*oG.ddxiMat;
    dfdxi.Ups1 = dia(3*oMXiSqDivGP2)*oG.ddxiMat;
    dfdxi.Ups2 = dia(-4*oMXiSqDivGP2)*oG.ddxiMat;
    dfdxi.UpsPlus = dia(2*xi.*divGP2);
    dfdxi.UpsMin = dia(-xi.*oMXiSqDivGP4)*oG.d2dxi2Mat ...
                   -dia(2*gDivP.*oMXiSqDivP2)*oG.d2dpdxiMat ...
                   -dia(2*xi.*gDivP3)*oG.ddpMat ...
                   +dia(2*divGP4 + 3*oMXiSqDivGP2)*oG.ddxiMat;            
    dfdxi.Pi = dia(-oMXiSqDivGP2)*oG.ddxiMat;

    %d2f/dpdxi term
    d2fdpdxi.UpsMin = dia(2*gamma.*oMXiSqDivP2)*oG.d2dpdxiMat ...
                      -dia(2*gDivP.*oMXiSqDivP2)*oG.ddxiMat;

    %f term (no derivative)
    f.Pi = dia(-gamma)*oG.d2dp2Mat ...
           -dia(divGP.*(2+3*p.*p))*oG.ddpMat ...
           -dia(oMXiSqDivGP2)*oG.d2dxi2Mat ...
           +dia(2*xi.*divGP2)*oG.ddxiMat;             

    %Save the precalculated parts   
    o.d2fdp2 = d2fdp2;
    o.dfdp = dfdp;
    o.d2fdxi2 = d2fdxi2;
    o.dfdxi = dfdxi;
    o.d2fdpdxi = d2fdpdxi;
    o.f = f;


    %%% Prefactors for the electron-ion operator
    o.preFactorEI = oN.Z*oN.nBar*oN.lnLambdaBar;
    o.ionD2fdxi2 = 0.5*gDivP.*oMXiSqDivP2;
    o.ionDfdxi   = -xi.*gDivP3;            
    o.ionDfdxi(end)   = 0;

end
