function UpdateKineticEquation(o,fOld,t,tOld)
    % Updates and builds the various coefficients needed for the
    % determination of the residual and Jacobian matrix in Newton's method.
    % These are used to solve for F and EHat consistently in the case
    % of an inductive E field.
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
    oN.heatSink.Update(fOld,t,tOld);

    %Assemble the matrix components
    oG = oN.grid;    
    dt = t-tOld;
    
    %Update the avalanche source and fetch
    %its result (if it exists)
    if ~isempty(oN.avalancheSource)
        oN.avalancheSource.Assemble(t, fOld);
        avaS = oN.avalancheSource.S * dt;
    else
        avaS = 0;
    end
    
    AHat = oN.sigma(t)*oN.synchrotronOperator...
           + oN.collisionOperator.C ...
           + oN.heatSink.sink;
    rhs = avaS;
    
    % Build and update particle source
    oN.particleSource.Update(t,tOld,fOld,AHat,rhs);
    
    o.inductiveCoefficients.R1 = oN.identityWithBoundaryConditions - dt*AHat;
    
    o.inductiveCoefficients.R2 = -dt*oN.EFieldOperator;
                %The minus sign is required, since EFieldOperator is
                %defined with a minus sign.
    
    %The old right-hand side is now in R3    
    o.inductiveCoefficients.R3 = -fOld - oN.particleSource.source - rhs;
    o.inductiveCoefficients.R3(end) = 0; %Neumann boundary condition at p=0
    o.inductiveCoefficients.R3(oG.idsPMax) = 0; %Dirichlet boundary condition at p=pMax
    
    o.inductiveCoefficients.R4 = oN.LHat(t)/dt * ...
                       ( oG.intdpdxi .* (oG.xiBig.*oG.pBig./oG.gammaBig) )';
    o.inductiveCoefficients.R5 = -oN.EHat(t)-o.inductiveCoefficients.R4*fOld;
            %EHat is the _applied_ external electric field!
end
