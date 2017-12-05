function Assemble(o,t) 
    % Assembles the precalculated collision operator parts into one
    % operator matrix. For the field-particle part of the electron-electron
    % operator, potentials calculated implicitly from the distribution are
    % used, so that the matrix can be applied repeatedly wihtout being
    % rebuilt. If the temperature or density have changed, the matrix does
    % however need to be rebuilt. t is the time.
    %
    % Usage:
    %   Assemble(t)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tBuild = tic;    

    oN     = o.norse;
    oG     = o.grid;
    th     = oN.Theta(t);
    invTh  = 1/th;
    invTh2 = invTh*invTh;
    invTh3 = invTh*invTh2;
    invTh4 = invTh*invTh3;
     
    p      = oG.pBig;    
    divP   = 1./p;    
    
    
    %Helper function for diagonal matrices, to reduce clutter
    dia = @(v) spdiags(v,0,oG.matSize,oG.matSize); 

    
    %%% Build the field-particle term %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if o.includeFPTerm
        %Get a general implicit operator for the potentials (this avoids
        %calculating the large final potential matrices before they are
        %actually needed)
        potMap = @(M) oN.legModesToBigVectorMap*M*oN.bigVectorToLegModeMap;

        %Calculate weighted matrices for the potentials 
        o.CalculateWeightedPotentialMatrices(t);

        %%% Add the various terms up as we go along to avoid saving a number of
        %%% large matrices in memory.
        % @@@ Is this efficient? @@@

        %Upsilon0 terms
        fPTerm = ( invTh*o.dUpsilon0.invTh + invTh2*o.dUpsilon0.invTh2 ) ...
                 * potMap(o.weightedU0Mat);

        %Upsilon2 terms
        fPTerm = fPTerm + ( invTh*o.dUpsilon2.invTh + invTh2*o.dUpsilon2.invTh2 )... 
                          * potMap(o.weightedU2Mat); 

        %UpsilonMin terms
        fPTerm = fPTerm...
               + (   invTh  * o.dUpsilonMin.invTh   ...
                   + invTh2 * o.dUpsilonMin.invTh2  ...
                   + invTh3 * o.dUpsilonMin.invTh3  ...
                   + invTh4 * o.dUpsilonMin.invTh4 )...
                 * potMap(o.weightedUMinMat); 

        %UpsilonPlus terms
        fPTerm = fPTerm + invTh2*o.dUpsilonPlus * potMap(o.weightedUPlusMat);

        %Pi terms
        fPTerm = fPTerm ...
               + (o.dPi.noTh + invTh*o.dPi.invTh)* potMap(o.weightedPiMat);

        %Finalize by including the prefactor
        fPTerm = o.preFactorFP(t) * fPTerm;
    else
        fPTerm = sparse(oG.matSize,oG.matSize);
    end
    
    %%% Build the test-particle term %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if o.includeTPTerm
        o.CalculateSpecialFunctions(t);

        tPTerm = dia(o.GPara)*oG.d2dp2Mat ...
               + dia( 2*divP.*o.GPara + o.dGParaDP + o.GK )*oG.ddpMat ...
               + dia(2*divP.*o.GK + o.dGKDP) ...               
               + dia(o.GPerp)*o.d2fdxi2 ...
               + dia(o.GPerp)*o.dfdxi;

       tPTerm  = o.preFactorTP(t) * tPTerm;
    else
       tPTerm  = sparse(oG.matSize,oG.matSize); 
    end
    
    %%% Build the electron-ion operator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ionD2Term   = dia(o.ionD2fdxi2)*oG.d2dxi2Mat;
    ionDTerm    = dia(o.ionDfdxi)*oG.ddxiMat;
    ionOperator = o.preFactorEI(t)*(ionD2Term + ionDTerm);

    
    %%% Build the combined operator (with boundary conditions) %%%%%%%%%%%%
    o.C = fPTerm + tPTerm + ionOperator + o.dampeningMat;
    
    idsLastRow = find(oG.ddpMat(end,:));
    o.C(end,:) = 0; %Get rid of some NaNs that shouldn't affect the result
    o.C(end,idsLastRow) = oG.ddpMat(end,idsLastRow); %Boundary condition at p=0    
    % Apply boundary condition at p=pMax. These lines are slow
    % compared to the other lines! 
    o.C(oG.idsPMax,:) = 0;                               
    o.C(sub2ind(size(o.C),oG.idsPMax,oG.idsPMax)) = 1; 
               % Accesses all elements corresponding to pMax on the
               % diagonal to enforce F(pMax)=0.      

    o.norse.timing.collOp = o.norse.timing.collOp + toc(tBuild);            
end
