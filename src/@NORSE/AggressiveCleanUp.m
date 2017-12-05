function AggressiveCleanUp(o)
    % Keep only the bare minimum which still makes the plots etc.
    % work, to reduce memory consumption. The object cannot be used
    % to continue the calculation!
    %
    % Usage: 
    %   AggresiveCleanUp()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    o.CleanUp();

    o.identityWithBoundaryConditions = [];            

    o.legModesToBigVectorMap = [];
    o.bigVectorToLegModeMap = [];
    o.timeAdvance.idsToSave = [];
    o.timeAdvance.timesToSave = [];
end
