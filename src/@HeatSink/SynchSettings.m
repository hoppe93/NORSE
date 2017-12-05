function SynchSettings(o)
    % Makes sure settings are consistent with those in NORSE.
    % Should be called before the calculation starts (or at restart).
    %
    % Usage:
    %   SynchSettings()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %We want these parameters inside the HeatSink object, since it
    %gives faster access, and we need them in every time step.

    %Flags
    oN = o.norse;            
    o.includeHeatSink = oN.includeHeatSink;
    o.includeTempChanges = ~oN.T.isScalar;
    o.calculateDensityRelatedChanges = ~oN.n.isScalar;
    o.includeDensityRelatedChanges = o.calculateDensityRelatedChanges && ...
                            oN.keepTempConstantWhenDensityChanges;
    o.enforceStrictHeatConservation = oN.enforceStrictHeatConservation;
    o.restrictHeatSinkRate = o.includeHeatSink*oN.maxHeatSinkRate;
    if oN.collisionOperatorMode == 1 && ~o.includeHeatSink
        o.includeHeatSink = true;
        oN.includeHeatSink = true; 
        warning('The heat sink has been automatically enabled, since it is required for the conservative linear collision operator.');
    end
    if ~o.includeHeatSink && ...
            (o.includeTempChanges || o.includeDensityRelatedChanges)
        o.includeHeatSink = true;
        oN.includeHeatSink = true;
        warning('The heat sink has been automatically enabled, since it is required for temperature or density changes.');
    end

    %Calculate an energy rate conversion factor            
    mc2 = oN.constants.m*oN.constants.c^2;
    o.rateNormalization = mc2*oN.fM0;            

    if o.restrictHeatSinkRate
        o.rateCutOff = o.restrictHeatSinkRate/o.rateNormalization;
        %This should be compared to k*SHInt
    end

    %Other
    o.preFactor = 1/(oN.ThetaBar*oN.kappaBar);            
end
