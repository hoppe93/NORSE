function Assemble(o,t,fOld)
    % Assembles the precalculated collision operator parts into one
    % operator matrix. Uses the potentials calculated from the
    % distribution at time t. This method also multiplies the source
    % operator S with the distribution function 'fOld' of the
    % previous timestep.
    %
    % Usage:
    %   Assemble(t)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    oN = o.norse;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the critical momentum, pc
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute Ec / Ec0 (for conversion to
    % instantaneous units)
    Ec = oN.CalculateEc(oN.T(t), oN.n(t)) / o.Ec0;
    EHat = oN.EHat(t);
    
    % Include inductive electric field if available
    if oN.timeAdvanceMode == 3
        EHatInd = oN.EHatInductive(oN.timeAdvance.iTimeStep-1);
        EHat = EHat + EHatInd;
    end
    EHat = EHat / Ec;
    
    pc = 1 / sqrt(EHat - 1);
    
    % Mask out the parts of the operator where p < pc,
    % which diverges and should not contribute to the
    % avalanche process anyway.
    
    o.S = (o.intOperator * fOld) / oN.lnLambdaBar(t);
    % Remove source at p < pc to avoid instabilities
    % due to the diverging nature of S.
    o.S(o.grid.pBig < pc) = 0;
    
end
