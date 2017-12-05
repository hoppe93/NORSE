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

    %Remove unused entries
    TrimSaveArrays@TimeAdvance(o,stepCounter,iSave); %Do all the things in the superclass 
    o.dtsUsed = o.dtsUsed(1:stepCounter-1); %Also modify dtsUsed
end
