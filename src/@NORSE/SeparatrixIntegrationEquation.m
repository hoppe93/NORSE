function dpdx = SeparatrixIntegrationEquation(xi,p,TNorm,EHat,...
                                              sigma,dPidp,dPidxi)
    % Describes the equation to integrate in the calculation of the
    % separatrix of the runaway region. The inputs are the
    % corresponding quantites in the method IntegrateSeparatrix.
    %
    % Usage:
    %   dpdx = SeparatrixIntegrationEquation(xi,p,TNorm,EHat,
    %                                           sigma,dPidp,dPidxi)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    g = sqrt(1+p*p);
    p2OMinXi2 = p/(1-xi*xi);
    numerator = EHat*xi + g*TNorm*ppual(dPidp,[p;xi]) ...
                                            - sigma*g*p*(1-xi*xi);
    denominator = EHat + TNorm/(g*p)*ppual(dPidxi,[p;xi]) ...
                                                    + sigma*xi*p/g;  
    dpdx = p2OMinXi2*numerator/denominator;            
end
