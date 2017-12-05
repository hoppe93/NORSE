function m = PrepareBoundaryConds(o,m,l)
    % Sets the proper l-dependent boundary conditions at p=0 and
    % p=pmax in the matrix m. l is the Legendre-mode index.
    %
    % Usage:
    %   m = PrepareBoundaryConds(m,l)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if l
        %l~=0 -- use Neumann condition at p=0                
        m(1,m(1,:)~=0) = 0; %This is faster for big sparse matrices (at least)
        m(1,1) = 1;                 
    else
        %l=0 -- use Dirichlet condition at p=0 
        nz1 = find(o.grid.ddp(1,:)); 
        m(1,nz1) = o.grid.ddp(1,nz1);
    end            
    m(end,m(end,:)~=0) = 0; 
    m(end,end) = 1;
end        
