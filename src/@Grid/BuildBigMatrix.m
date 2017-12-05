function m = BuildBigMatrix(o,fOfXi,fOfP,doP0Dependence)
    % Construct a matrix that can be used for operating on a
    % function on the grid (in its vector representation). Gives a
    % vector when multiplying a function on the grid.
    %
    % Usage:
    %   m = BuildBigMatrix(fOfXi,fOfP,doP0Dependence)
    %
    % fOfXi describes the xi dependence and is a matrix with
    % dimenstion [nXi,nXi]. Similarly, fOfP describes the p
    % dependence and is [nP,nP]. The final argument determines
    % whether to perform the operations necessary to get the
    % correct beahvior at the special point p=0. If there is no
    % p dependence (i.e. fOfP = speye(nP,nP)), 0 can be passed here.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    m = sparse(o.matSize,o.matSize); 
    m(1:(end-1),1:(end-1)) = kron(fOfXi,fOfP(2:end,2:end)); 
                                    %Everything not related to p=0
    if doP0Dependence
        m(1:(end-1),end) = kron(ones(o.nXi,1),fOfP(2:end,1)); 
                            %Other grid points' depedence on p=0  
        xi0Ids = (o.xi0Id-1)*(o.nP-1)+(1:(o.nP-1));    
        m(end,xi0Ids) = fOfP(1,2:end); 
                %p=0's dependence on other grid points (at xi=0, 
                %where we specify the boundary condition)    
        m(end,end) = fOfP(1,1); %p=0's dependence on p=0  
    end
end
