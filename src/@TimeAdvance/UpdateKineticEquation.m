function UpdateKineticEquation(o,fOld,t,tOld)
    % Updates and builds the matrix and right-hand side describing
    % the kinetic equation.
    %
    % Usage:
    %   UpdateKineticEquation(fOld,t,tOld)
    %
    % fOld and tOld are the distribution and time at the previous
    % time step. t is the time at the current time step.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    oN = o.norse;
    
    %Calculate the potentials from f and update the collision operator 
    fls  = oN.MapBigVectorToLegModes(fOld);
    oN.potentials.Update(fls); 
    oN.collisionOperator.Assemble(t);

    %Update the heat and particle sources
    oN.particleSource.Update(t,tOld);
    oN.heatSink.Update(fOld,t,tOld);

    %Assemble the matrix (including time-advancement scheme)                
    dt = t-tOld;            
    o.matrix = oN.identityWithBoundaryConditions - ...
               dt*(oN.EHat(t)*oN.EFieldOperator ...
                   + oN.sigma(t)*oN.synchrotronOperator...
                   + oN.collisionOperator.C ...
                   + oN.heatSink.sink); %This is about 5 times faster 
                                       %than calling sparse() on
                                       %index vectors.

    %Build the right-hand side (including the particle source)
    o.rhs = fOld + oN.particleSource.source;
    o.rhs(end) = 0; %Neumann boundary condition at p=0
    o.rhs(oN.grid.idsPMax) = 0; %Dirichlet boundary condition at p=pMax
end
