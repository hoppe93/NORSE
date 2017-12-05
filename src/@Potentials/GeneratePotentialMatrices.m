function GeneratePotentialMatrices(o)
    % Generates the matrices needed to efficiently calculate the
    % potentials from the distribution.
    %
    % Usage:
    %   GeneratePotentialMatrices()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    p = o.grid.p;
    nP = o.grid.nP;
    gamma = o.grid.gamma;
    gamma2 = gamma.*gamma;
    nL = o.norse.nL;

    %Helper function for diagonal matrices, to reduce clutter
    dia = @(v) spdiags(v,0,nP,nP);
    o.La_lBase = dia(gamma2)*o.grid.d2dp2 + dia(2./p+3*p)*o.grid.ddp;

    o.CalculateJAndY();            

    %Get the mapping matrices from f to the various potentials for
    %each Legendre mode l
    sPotentials(nL) = o.GeneratePotentialMatricesForL(nL-1);
    for l = 0:(nL-2)                
        sPotentials(l+1) = o.GeneratePotentialMatricesForL(l);
    end

    %Assemble into one big integrator matrix for each potential,
    %which we can apply to a vector with structure
    % [f_0(p_1,...,p_nP);f_1(p_1,...,p_nP);...;f_nL(p_1,...,p_nP),f_0(p0)].
    %The matrices in sPotential have size (nP-1)x(nP-1) (the first
    %grid point is excluded). On the last row we enforce that
    %dPsi_0/dp=0 at p=0 by explicitly determining Psi_0(p=0) from
    %the expression for the finite difference derivative and the
    %(known) value of Psi_0 at grid points p=1,...,4.
    cPotentials = struct2cell(sPotentials); 
                    %This gives a cell array with has the fields 
                    %Upsilon0,...,Pi1 for a given l in each column 
    ddp = o.grid.ddp(1,:);            
    nStencil = nnz(ddp);  
    ddpBC = -ddp(2:nStencil)/ddp(1);

    o.Upsilon0Matrix = blkdiag(cPotentials{1,:},0);
    o.Upsilon0Matrix(end,:) = ddpBC*o.Upsilon0Matrix(1:(nStencil-1),:);
    %For a 5-point stencil, the latter is equivalent to:
    %o.Upsilon0Matrix(end,:) = -1/ddp(1)*( ddp(2)*o.Upsilon0Matrix(1,:)...
    %                                     +ddp(3)*o.Upsilon0Matrix(2,:)...
    %                                     +ddp(4)*o.Upsilon0Matrix(3,:)...
    %                                     +ddp(5)*o.Upsilon0Matrix(4,:) ); 
    % where o.PsiMatrix(i,:)*F ==> Psi(p(i+1)), for a potential Psi.

    o.Upsilon1Matrix = blkdiag(cPotentials{2,:},0);
    o.Upsilon1Matrix(end,:) = ddpBC*o.Upsilon1Matrix(1:(nStencil-1),:);

    o.Upsilon2Matrix = blkdiag(cPotentials{3,:},0);
    o.Upsilon2Matrix(end,:) = ddpBC*o.Upsilon2Matrix(1:(nStencil-1),:);

    Pi0Matrix = blkdiag(cPotentials{4,:},0);
    Pi0Matrix(end,:) = ddpBC*Pi0Matrix(1:(nStencil-1),:);

    Pi1Matrix = blkdiag(cPotentials{5,:},0);
    Pi1Matrix(end,:) = ddpBC*Pi1Matrix(1:(nStencil-1),:);

    o.PiMatrix = 2*Pi1Matrix - Pi0Matrix; 
end        
