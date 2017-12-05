function BuildAndCalculateMagnitude(o,t,tOld)
    % Builds the particle source operator and calculates the
    % appropriate magnitude to perform specified changes to the
    % electron density.
    %
    % Usage:
    %   BuildAndCalculateMagnitude(t,tOld)
    %
    % t is the time at the current time step and tOld is the time
    % at the previous time step.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~o.norse.n.isScalar
        %Check for density changes
        oN = o.norse;
        dn = oN.n(t)-oN.n(tOld);
        if dn
            %The density has changed, calculate the needed source magnitude
            o.densityChangeMagnitude = dn/o.particleSourceDensityMoment(t);

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
end
