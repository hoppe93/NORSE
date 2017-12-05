function TrimSaveArrays(o,stepCounter,iSave)
    % Removes unused space at the end of the vectors/matrices used 
    % for saving time step data. This is used with the adaptive
    % time-step scheme or when aborting the calculation prematurely.
    % 
    % Usage: 
    %   TrimSaveArrays(stepCounter,iSave)
    %
    % stepCounter is the total number of steps taken so far and
    % iSave is the save step index (the current position in the 
    % save arrays).
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    oN = o.norse;
    
    %Remove unused entries
    o.gmresFlags = o.gmresFlags(1:stepCounter-1);
    o.gmresRess  = o.gmresRess(1:stepCounter-1);
    o.gmresIters = o.gmresIters(1:stepCounter-1);            
    o.allTimes   = o.allTimes(1:stepCounter-1);
    o.nTimeSteps = numel(o.allTimes);    

    o.times                 = o.times(1:iSave-1);
    oN.f(:,iSave:end)        = [];
    if oN.savePotentials
        oN.Upsilon0(:,iSave:end) = [];
        oN.Upsilon1(:,iSave:end) = [];
        oN.Upsilon2(:,iSave:end) = [];
        oN.Pi(:,iSave:end)       = [];
    end
    oN.energyChangeMagnitudes(iSave:end) = [];
    o.idsToSave(iSave:end)  = [];
    if ~isempty(oN.EHatInductive)
        oN.EHatInductive(iSave:end)      = []; 
    end
    o.nSaveSteps            = iSave-1; %Make sure everything is consistent
end
