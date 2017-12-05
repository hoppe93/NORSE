function s = GenerateBoundaryConditionsForInt(o,l)
    % Uses Eq. (31) in Braams & Karney [PoF B, 1 (1989), 1355] to
    % calculate boundary conditions for the differential
    % calculation of the potentials at p=pMax. We need the quantity
    % N*p^2/gamma in the equation for the boundary condition.
    %
    % Usage:
    %   s = GenerateBoundaryConditionsForInt(l)
    %
    % l is the Legendre-mode index and s is a struct with the
    % fields: Upsilon0, Upsilon1, Upsilon2, Pi0 and Pi1.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    p2g = o.grid.p.^2./o.grid.gamma;
    NL0 = o.GetNLABoundary(l,'0');
    NL02 = o.GetNLABoundary(l,'02');
    NL022 = o.GetNLABoundary(l,'022');
    NL1 = o.GetNLABoundary(l,'1');
    NL11 = o.GetNLABoundary(l,'11');

    %At p=pMax, the second integral vanishes, so that the integrand
    %is just N*p^2/gamma
    s = struct();
    s.Upsilon0 = p2g'.*NL0; 
    s.Upsilon1 = p2g'.*NL02;
    s.Upsilon2 = p2g'.*NL022;
    s.Pi0 = p2g'.*NL1;
    s.Pi1 = p2g'.*NL11; 
end
