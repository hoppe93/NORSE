function InitializeParts(o) 
    % Calculates the time-independent prefactors and grid functions that
    % multiply the various terms in the collision operator.
    %
    % Usage:
    %   InitializeParts()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    o.SynchSettings();
    
    %%% Prefactors for the electron-electron operator
    if o.includeTPTerm
        o.InitializeTestParticleTerm();
    end
    if o.includeFPTerm
        o.InitializeFieldParticleTerm();
    end
    
    %%% Prefactors for the electron-ion operator
    oN = o.norse;
    oG = o.grid;
    p     = oG.pBig;
    gamma = oG.gammaBig;
    xi    = oG.xiBig;
    
    divP   = 1./p;
    divP2  = divP.*divP;
    gDivP  = gamma.*divP;
    gDivP3 = gDivP.*divP2;    
    oMXiSqDivP2 = (1-xi.*xi).*divP2;
        
    o.preFactorEI = oN.Z*oN.nBar*oN.lnLambdaBar;
    o.ionD2fdxi2 = 0.5*gDivP.*oMXiSqDivP2;
    o.ionDfdxi   = -xi.*gDivP3;            
    o.ionDfdxi(end)   = 0;
end
