function iSave = CheckAbortCriterion(o,fOld,t,isSaveStep,iSave,iteration)
    % Check whether to abort the calculation because the system has
    % entered the slide-away regime.
    %
    % Usage:
    %   iSave = CheckAbortCriterion(fOld,t,isSaveStep,iSave,iteration)
    %
    % iSave is the save step index, fOld is the distribution at the
    % previous time step and t is the time at the current time
    % step. isSaveStep is a boolean, and iteration is the time step
    % index.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    doAbortS = 0;
    doAbortG = 0;
    
    oN = o.norse;
    
    %Check whether the separatrix has dissappeared                    
    if oN.abortWhenOnlyRunaways && ~oN.FindParallelPCrit('inst',t)
        doAbortS = 1;
        oN.Print('\n         Aborting the run due to a transition to the slide-away regime!');
    end    
    if ~isempty(oN.GeneralAbortFunction)
        %Call the user-supplied function to determine if an abort criterion
        %has been fulfilled
        doAbortG = oN.GeneralAbortFunction(oN,t,isSaveStep,iSave,iteration);
        if doAbortG
            oN.Print('\n         Aborting the run since a user-specified criterion is fulfilled!');
        end
    end    
    if doAbortS || doAbortG
       %Abort. Make sure the final distribution is saved.
        if ~isSaveStep
            o.SaveStepData(iSave,fOld,t);                                        
            iSave = iSave+1;
        end
        o.TrimSaveArrays(iteration+1,iSave);
        o.abortFlag = true;            
    end
end
