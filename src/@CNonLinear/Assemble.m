function Assemble(o,t) 
    % Assembles the precalculated collision operator parts into one
    % operator matrix. Uses the potentials calculated from the
    % distribution at time t.
    %
    % Usage:
    %   Assemble(t)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    tBuild = tic;    

    %Helper function for diagonal matrices, to reduce clutter
    dia = @(v) spdiags(v,0,o.grid.matSize,o.grid.matSize); 

    %Get the potentials
    Ups0 = o.potentials.Upsilon0;
    Ups1 = o.potentials.Upsilon1;
    Ups2 = o.potentials.Upsilon2;            
    UpsPlus = o.potentials.UpsilonPlus;
    UpsMin = o.potentials.UpsilonMinus;
    Pi = o.potentials.Pi;            

    %%% Each derivative of f has a prefactor that depends only on
    %%% the potentials. Calculate that prefactor and attach the
    %%% appropriate differentiation matrix (to apply to f).

    %%% Electron-electron collision operator
    d2fdp2Term   = dia( o.d2fdp2.Ups0*Ups0 ...
                       +o.d2fdp2.Ups2*Ups2 ...
                       +o.d2fdp2.UpsMin*UpsMin)...
                   *o.grid.d2dp2Mat;

    dfdpTerm     = dia( o.dfdp.Ups0*Ups0 ... 
                       +o.dfdp.Ups1*Ups1 ...
                       +o.dfdp.Ups2*Ups2 ...
                       +o.dfdp.UpsMin*UpsMin...
                       +o.dfdp.Pi*Pi)...
                   *o.grid.ddpMat;                   

    d2fdxi2Term  = dia( o.d2fdxi2.UpsPlus*UpsPlus...
                       +o.d2fdxi2.UpsMin*UpsMin)...
                   *o.grid.d2dxi2Mat;

    dfdxiTerm    = dia( o.dfdxi.Ups0*Ups0 ...
                       +o.dfdxi.Ups1*Ups1 ...
                       +o.dfdxi.Ups2*Ups2 ...
                       +o.dfdxi.UpsPlus*UpsPlus...
                       +o.dfdxi.UpsMin*UpsMin...
                       +o.dfdxi.Pi*Pi)...
                   *o.grid.ddxiMat; 

    d2fdpdxiTerm = dia(o.d2fdpdxi.UpsMin*UpsMin)...
                   *o.grid.d2dpdxiMat;

    fTerm        = dia(o.f.Pi*Pi);


    %%% Electron-ion collision operator            
    ionD2Term = dia(o.ionD2fdxi2)*o.grid.d2dxi2Mat;
    ionDTerm = dia(o.ionDfdxi)*o.grid.ddxiMat;

    %%% Build the combined operator                
    C = o.preFactorEE(t)*(  d2fdp2Term + dfdpTerm + d2fdxi2Term ...
                          + dfdxiTerm + d2fdpdxiTerm + fTerm   )...
        + o.preFactorEI(t)*(ionD2Term + ionDTerm) ...
        + o.dampeningMat;
    idsLastRow = find(o.grid.ddpMat(end,:));            
    C(end,idsLastRow) = o.grid.ddpMat(end,idsLastRow); %Boundary condition at p=0    
    % Apply boundary condition at p=pMax. These lines are slow
    % compared to the other lines! 
    C(o.grid.idsPMax,:) = 0;                               
    C(sub2ind(size(C),o.grid.idsPMax,o.grid.idsPMax)) = 1; 
               % Accesses all elements corresponding to pMax on the
               % diagonal to enforce F(pMax)=0.
    o.C = C;             

    o.norse.timing.collOp = o.norse.timing.collOp + toc(tBuild);            
end
