function CalculateNumberOfSteps(o)
    %Perform some initialization needed to set up the rest of the NORSE
    %calculation. Most of the actual time-advance-specific initialization
    %is done in InitializeTimeStepArrays()
        
    oN = o.norse;
    
    if isempty(oN.nTimeSteps)
        oN.nTimeSteps = ceil((oN.tMax-eps)/oN.dt)+1;
        if ~oN.nSaveSteps                 
            %Save all steps
            oN.nSaveSteps = oN.nTimeSteps; 
        else
            %Save some steps (unless fewer steps than specified
            %are actually taken)
            oN.nSaveSteps = min(oN.nSaveSteps,oN.nTimeSteps); 
        end
        switch oN.timeStepSaveMode
            case 0 %linear
                o.idsToSave = unique(round(linspace(1,oN.nTimeSteps,oN.nSaveSteps)));
            case 1 %logarithmically increasing
                o.idsToSave = unique(round(logspace(0,log10(oN.nTimeSteps),oN.nSaveSteps)));
            otherwise
                error('Invalid timeStepSaveMode.');
        end
        %Some may have been removed by unique()
        oN.nSaveSteps = numel(o.idsToSave);
    else
        %Restart, nothing to do here. It will be taken care of in
        %InitializeTimeStepArrays()
    end
end