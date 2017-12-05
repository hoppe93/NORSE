function sPot = GeneratePotentialMatricesForL(o,l)
    % Calculates the matrices needed to calculate the various
    % potentials from the distribution for a given Legendre mode l.
    % The potentials are returned in a struct sPot with fields 
    %   Upsilon0, Upsilon1, Upsilon2, Pi0 and Pi1.
    %
    % Usage:
    %   sPot = GeneratePotentialMatricesForL(l)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    nP = o.grid.nP;
    p2 = o.grid.p.*o.grid.p;
    invP2 = 1./(o.grid.p.*o.grid.p);

    %Define anonynoums functions to take care of the a-dependence
    %of the differential operator L_a
    if l 
        La = @(a) o.PrepareBoundaryConds(o.La_lBase +... 
                          spdiags(1-a^2-l*(l+1).*invP2,0,nP,nP), l);
    else                
        La = @(a) o.PrepareBoundaryConds(...
                            o.La_lBase +  (1-a^2)*speye(nP,nP), l);
    end

    %Calculate boundary conditions for the potentials at p=pmax           
    sBounds = o.GenerateBoundaryConditionsForInt(l); 


    %%% Solve for the operator matrices for each potential
    %%% according to the differential relations between them

    IBase = eye(nP); %No point in making this sparse, since the 
                     %potentials matrices will be full anyway
    if l
        IBase(1,1) = 0;            
    end

    %Upsilon0_l
    rhs = IBase;
    rhs(end,:) = o.grid.pWeights'.*sBounds.Upsilon0; %Boundary condition
    Upsilon0 = La(0)\rhs; %This is essentially the inverse of La(0).
                          %We can multiply it with f_l to find
                          %Upsilon_0 for aritrary f.

    %Upsilon1_l
    rhs = Upsilon0;
    rhs(1,:) = 0;
    rhs(end,:) = o.grid.pWeights'.*sBounds.Upsilon1; 
    Upsilon1 = La(2)\rhs;    

    %Upsilon2_l
    rhs = Upsilon1;
    rhs(1,:) = 0;
    rhs(end,:) = o.grid.pWeights'.*sBounds.Upsilon2; 
    Upsilon2 = La(2)\rhs;

    %Pi0_l
    rhs = IBase;
    rhs(1,1) = 0;
    rhs(end,:) = o.grid.pWeights'.*sBounds.Pi0; 
    Pi0 = La(1)\rhs;

    %Pi1_l
    rhs = Pi0;
    rhs(1,:) = 0;
    rhs(end,:) = o.grid.pWeights'.*sBounds.Pi1; 
    Pi1 = La(1)\rhs;            


    %Do not return the first row. The behavior there is known (the
    %potential=0 for l>0, and the derivative of the potential=0 for
    %l=0), and this structure makes it easier to construct an
    %operator to act on the entire f at once. Similarly, do not
    %return the first column.
    sPot = struct('Upsilon0',sparse(Upsilon0(2:end,2:end)),...
                  'Upsilon1',sparse(Upsilon1(2:end,2:end)),...
                  'Upsilon2',sparse(Upsilon2(2:end,2:end)),...
                  'Pi0',sparse(Pi0(2:end,2:end)),...
                  'Pi1',sparse(Pi1(2:end,2:end)));   
end
