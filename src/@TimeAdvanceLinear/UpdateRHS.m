function UpdateRHS(o,fOld,t,tOld)
    % Updates the right-hand side of the kinetic equation, which
    % changes in every time step.
    %
    % Usage:
    %   UpdateRHS(fOld,t,tOld)
    %
    % fOld and tOld are the distribution and time at the previous
    % time step. t is the time at the current time step.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    oN = o.norse;
    
    %Update the heat and particle sources
    oN.particleSource.Update(t,tOld);
    oN.heatSink.Update(fOld,t,tOld);
    
    %Build the right-hand side (including the particle source)
    o.rhs = fOld + oN.particleSource.source - oN.heatSink.sink*fOld;
    o.rhs(end) = 0; %Neumann boundary condition at p=0
    o.rhs(oN.grid.idsPMax) = 0; %Dirichlet boundary condition at p=pMax
end
