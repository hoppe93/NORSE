function InitializeTestParticleTerm(o) 
    % Calculates the time-independent prefactors for the test-particle term
    % of the electron-electron collision operator.
    %
    % Usage:
    %   InitializeTestParticleTerm()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Preliminaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    oN = o.norse;
    oG = o.grid;
    p     = oG.pBig;    
    xi    = oG.xiBig;

    %Helper function for diagonal matrices, to reduce clutter
    dia = @(v) spdiags(v,0,oG.matSize,oG.matSize); 
    
    %Miscellaneous grid quantities
    divP  = 1./p;
    divP2 = divP.*divP;
    oMXiSq  = 1-xi.*xi;
    oMXiSqDivP2 = oMXiSq.*divP2;
    
    %Prefactor
    o.preFactorTP = oN.nBar*oN.lnLambdaBar;
    
    %xi-derivative terms
    o.d2fdxi2 =  dia(oMXiSqDivP2)*oG.d2dxi2Mat;
    o.dfdxi   = -dia(2*xi.*divP2)*oG.ddxiMat;    
end
