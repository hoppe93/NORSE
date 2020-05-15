function BuildAndCalculateMagnitude(o,t,tOld,fOld,matrix,rhs)
    % Builds the particle source operator and calculates the
    % appropriate magnitude to perform specified changes to the
    % electron density.
    %
    % Usage:
    %   BuildAndCalculateMagnitude(t,tOld,fOld,matrix,rhs)
    %
    % t is the time at the current time step and tOld is the time
    % at the previous time step.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    oN = o.norse;
    
    if oN.conservativeParticleSource == 0
        if ~oN.n.isScalar
            %Check for density changes
            dn = oN.n(t)-oN.n(tOld);
            if dn
                %The density has changed, calculate the needed source magnitude
                psDensityMoment = o.particleSourceDensityMoment(t);
                o.densityChangeMagnitude = dn/psDensityMoment;

                %Build the operator                    
                fM = oN.maxwellianPreFactor(t)*exp(-o.gammaMin1/oN.Theta(t));                
                o.particleSourceOperator = ...
                      ( o.gammaMin1/oN.Theta(t) + o.particleSourceAP(t) ).*fM;
                        %This dependes on density (and temperature),
                        %and so is always different when it's needed!
                        %No point in pre-calculating.
            else
                o.densityChangeMagnitude = 0;
            end
        else
            o.densityChangeMagnitude = 0;
        end
    % Approximately conservative source
    elseif oN.conservativeParticleSource == 1
        nOld = o.grid.intdpdxi' * (oN.fM0 * fOld);
        dn = oN.n(t) - nOld;
        dt = t-tOld;
        psDensityMoment = o.particleSourceDensityMoment(t);
        
        % Calculate the density moment of the kinetic equation terms
        kineqDensityMoment = oN.fM0 * o.grid.intdpdxi' * (matrix * fOld + rhs);
        
        %Build the operator                    
        fM = oN.maxwellianPreFactor(t)*exp(-o.gammaMin1/oN.Theta(t));                
        o.particleSourceOperator = ...
              ( o.gammaMin1/oN.Theta(t) + o.particleSourceAP(t) ).*fM;
        %This dependes on density (and temperature),
        %and so is always different when it's needed!
        %No point in pre-calculating.
        
        % Update the density change magnitude
        o.densityChangeMagnitude = (dn - dt*kineqDensityMoment) / psDensityMoment;
    else
        error('Invalid particle source specified.');
    end
end

