function iSave = SaveState(o,iSave,fOld,t)
    % Saves the state of the calculation and optionally plots the
    % parallel distribution.
    %
    % Usage:
    %   iSave = SaveState(iSave,fOld,t)
    %
    % iSave is the index into the save arrays (which is modified).
    % fOld is the distribution at the previous time step and t is
    % the time at the current time step.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    o.SaveStepData(iSave,fOld,t);                                        
    iSave = iSave+1;
    %Update "live" plot
    if o.norse.show1DTimeEvolution
        o.norse.Refresh1DTimeEvolutionPlot(fOld)
    end            
end
