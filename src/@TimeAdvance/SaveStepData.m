function SaveStepData(o,iSave,f,t)
    % Saves the state of the solution at a given time step for
    % post-processing and visualization.
    %
    % Usage:
    %   SaveStepData(iSave,f,t)
    %
    % iSave is the save step index (the position in the save
    % arrays), f is the distribution and t is the time.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    oN = o.norse;
    
    oN.f(:,iSave) = f;
    if oN.savePotentials
        oN.Upsilon0(:,iSave) = oN.potentials.Upsilon0;
        oN.Upsilon1(:,iSave) = oN.potentials.Upsilon1;
        oN.Upsilon2(:,iSave) = oN.potentials.Upsilon2;
        oN.Pi(:,iSave) = oN.potentials.Pi;
    end
    oN.energyChangeMagnitudes(iSave) = oN.heatSink.energyChangeMagnitude;
    oN.densityChangeMagnitudes(iSave) = ...
                            oN.particleSource.densityChangeMagnitude;
    o.times(iSave) = t;
end
