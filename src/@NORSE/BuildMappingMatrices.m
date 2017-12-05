function BuildMappingMatrices(o)
    % Constructs matrices for mapping between the 2D
    % finite-difference grid and finite-difference--Legendre-mode
    % representations.
    %
    % Usage:
    %   BuildMappingMatrices()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nP  = o.grid.nP;
    nXi = o.grid.nXi;

    %%%Build a matrix describing the mapping from a Legendre-mode
    %%%vector to a vector on the 2D finite-difference grid
    mapMat = sparse((nP-1)*nXi+1, (nP-1)*o.nL+1); 
    cellOfRowBlocks = cell(nXi,1);

    %Make a (nP-1)x(nP) block with ones on the first upper
    %diagonal, which we will then duplicate horizontally, weighted
    %by P_l, to form a "row block". Row block i gives the mapping
    %to f(p,xi_i).
    identityBlock = speye(nP-1);                        

    Pls = o.LegendrePolynomials(o.nL-1,o.grid.xi');            
    for iXi = 1:nXi
        PlsAtXi = Pls(:,iXi)';
        cellOfRowBlocks{iXi} = kron(PlsAtXi,identityBlock);
    end

    mapMat(1:(end-1),1:(end-1)) = cell2mat(cellOfRowBlocks);            
    mapMat(end,end) = 1; %Handle f(p=0)                  
    o.legModesToBigVectorMap = mapMat;            


    %%%Invert the matrix to get the inverse mapping from the
    %%%2D finite-difference grid
    %%%vector to Legendre modes
    %This should give accuracy to machine precision. The matrix is
    %rectangular, so we need to use the pseudo inverse. The
    %following is equivalent to pinv(full(mapMat)), but works for
    %sparse matrices (and the result is pretty sparse, so storing
    %is not an issue):
    r = qr(mapMat,0);
    o.bigVectorToLegModeMap = r\(r'\mapMat');            
end 
