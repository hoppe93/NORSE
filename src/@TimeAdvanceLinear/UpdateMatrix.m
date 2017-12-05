function UpdateMatrix(o,fOld,t,tOld)
    % Rebuilds the matrix describing the E-field, synchrotron and collision
    % operator terms. Only needs to be called when parameters change.    
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
    

    %Assemble the matrix (including time-advancement scheme)                
    dt = t-tOld;            
    o.matrix = oN.identityWithBoundaryConditions - ...
               dt*(oN.EHat(t)*oN.EFieldOperator ...
                   + oN.sigma(t)*oN.synchrotronOperator...
                   + oN.collisionOperator.C); %This is about 5 times faster 
                                       %than calling sparse() on
                                       %index vectors.
end
