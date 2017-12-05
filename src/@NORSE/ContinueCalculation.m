function ContinueCalculation(o)
    % Restart/continue the calculation. Use existing settings
    % and/or time-dependent parameters (the user may change
    % settings before this call). tMax signifies to
    % _total_ calculation time, not the time to run after
    % restarting. Similarly, nSaveSteps is for the entire run
    % (although no previously saved steps will be removed, and the
    % final time step will always be saved). When restarting, the
    % timeStepSaveMode setting will be ignored and linear spacing
    % will be used.
    %
    % Timing in the final print out may be erroneous if the
    % time-advancement scheme is changed before restart.
    %
    % !!! Disclaimer: The user is responsible for making sure that
    % updated time-dependent parameters remain consistent with
    % already saved time steps! Similarly, things that change the
    % normalization (such as modifying o.referenceTime) invalidate
    % the entire run. !!!
    %
    % Usage:
    %   ContinueCalculation()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Make sure the grid hasn't changed            
    if size(o.f,1) ~= ( (o.grid.nP-1)*o.grid.nXi +1 )
        error('Grid settings have been changed -- this is not allowed.'); 
    end
    if o.tMax == o.times(end)
        warning('tMax has already been reached. No calculation will be performed.');
        return
    end
    nOld = numel(o.times);
    
    %Check for changes to savePotentials
    if o.savePotentials && isempty(o.Pi)
        %User has enabled savePotentials -- issue a warning
        warning('Cannot save the potentials. They must be saved throughout the entire run.');
        o.savePotentials = 0;
    elseif ~o.savePotentials && ~isempty(o.Pi)
        %User has disabled savePotentials -- clear the fields
        warning('Removing previously saved potentials.');
        o.Upsilon0 = [];
        o.Upsilon1 = [];
        o.Upsilon2 = [];
        o.Pi       = [];
    end
    if o.runawayRegionMode==3 && ~o.savePotentials
	error('runawayRegionMode=3 requires that the potentials be saved.')
    end

    o.ProcessPhysicalParameters(); %In case parameters have changed or been updated
    if isempty(o.potentials)
        %Rebuild things that have been cleared by CleanUp().
        o.RebuildInternalVariables();                
    end
    
    %Generate new time-advance and collision operator objects to allow for
    %changing time-advance methods and collision operators. The
    %time-advance object will inherit all common properties from the
    %existing object.
    o.InitializeTAandC();
    o.collisionOperator.InitializeParts();

    tStart = tic;
    o.timeAdvance.AdvanceInTime(true); %Run            
    o.PostProcess(nOld+1);            
    o.timing.total = o.timing.total + toc(tStart);
    o.PrintTimings();
end
