function InitializeArtificialDamping(o) 
    % Perpare an artifical dampening term for the boundary condition at
    % p_max. Can be used to reduce ringing in the tail of the distribution.
    % The effect is controlled by the NORSE properties dampeningStrength
    % and dampeningWidth.
    %
    % Usage: 
    %   InitializeArtificialDamping()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    oN = o.norse;
    oG = o.grid;
    width = oN.dampeningWidth*oG.pMax; %dampeningWidth in per cent
    fakeViscosity = exp((oG.p-oG.pMax)/width);
    dampeningVector = oN.dampeningStrength*fakeViscosity;
    o.dampeningMat = oG.BuildBigMatrix(speye(oG.nXi),...
                        spdiags(dampeningVector,0,oG.nP,oG.nP)...
                        *oG.d2dp2,0);
end
