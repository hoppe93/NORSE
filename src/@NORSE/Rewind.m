function Rewind(o,iSave)
    % Make an intermediate save step iSave the final step, i.e. "rewind"
    % the calculation to a previous state. Can be used together with
    % ContinueCalculation() to (for instance) redo part of a run with a
    % smaller time step.
    %
    % Usage:
    %   Rewind(iSave)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Check that the save step is legit
    if iSave > o.nSaveSteps
        fprintf('Invalid time step. Cannot go back to the future.\n');
        return
    elseif iSave == o.nSaveSteps
        fprintf('No need to rewind -- the calculation is already in the desired state.\n');
        return
    elseif iSave <= 1
        fprintf('Rewinding to the first time step does not make sense. Do a new run instead.\n');
        return
    else
        %Nothing to do
    end
    
    %Get all property names
    oTA = o.timeAdvance;
    props  = properties(o);    
    propsTA = properties(oTA);
    nSave  = o.nSaveSteps;
    nSteps = o.nTimeSteps;
    iStep  = find(oTA.allTimes == o.times(iSave),1);
    if isempty(iStep)
        error('The NORSE object seems to be corrupted.');
    end
    
    %Loop through the relevant properties and remove the data to discard
    for i = 1:numel(props)
        prop       = props{i};
        propLength = size(o.(prop),2);
        if propLength == nSave
            o.(prop)(:,(iSave+1):end) = [];
        end
    end
    
    %Do the same for the timeAdvance object
    for i = 1:numel(propsTA)
        prop       = propsTA{i};
        propLength = size(oTA.(prop),2);
        if propLength == nSave
            oTA.(prop)(:,(iSave+1):end) = [];
        elseif propLength == nSteps
            oTA.(prop)(:,(iStep+1):end) = [];
        end
    end
    
    %Make sure the settings are consistent    
    oTA.nSaveSteps = iSave;    
    oTA.nTimeSteps = iStep;    
    oTA.SynchAfterCalculation();
    
    o.tMax   = oTA.allTimes(end);
    oTA.tMax = oTA.allTimes(end);
    
    oTA.abortFlag = 0; %Reset the abort flag in case it has been triggered 
    
    o.ResetTimings();  %Clear the timers since they will be incorrect    
    o.Print('Successfully rewinded the calculation to time step %d.\n',iSave);
end
