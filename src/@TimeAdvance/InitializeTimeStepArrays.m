function [fOld,iSave,firstStep] = InitializeTimeStepArrays(o,varargin)
    %Initializes the time-step vector and various save arrays for
    %the deterministic time-advancement schemes. Also handles restarts.
    %
    % Usage:
    %    [fOld,iSave,firstStep] = InitializeTimeStepArrays()
    %    [fOld,iSave,firstStep] = InitializeTimeStepArrays(true) -- restart
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    oN = o.norse;
    if nargin == 2 && varargin{1}
        %Restart -- extend arrays rather than initializing them

        %Extract information from the existing save arrays
        fOld            = oN.f(:,end);                        
        nOld            = o.nTimeSteps;                        
        firstStep       = nOld+1;
        nOldSave        = numel(o.times);
        iSave           = nOldSave + 1;
        tOld            = o.times(end);

        %Add new time steps
        nNew = ceil((o.tMax-tOld-eps)/oN.dt);
        o.nTimeSteps = nOld + nNew;
        times    = linspace(tOld,o.tMax,nNew + 1); 
                    %The +1 is to get the correct spacing;
                    %we will not use the first point
        o.allTimes   = [o.allTimes,times(2:end)];

        o.gmresFlags    = [o.gmresFlags,zeros(1,nNew)];
        o.gmresRess     = [o.gmresRess,zeros(1,nNew)];
        o.gmresIters    = [o.gmresIters,zeros(1,nNew)];                        

        %Determine how many and what steps to save
        if o.nSaveSteps > o.nTimeSteps
            %Save all
            ids = (o.idsToSave(end)+1):o.nTimeSteps;
            o.idsToSave = [o.idsToSave,ids];            
        elseif o.nSaveSteps <= nOldSave
            %Save only the final step (and keep all the old
            %time steps)
            o.idsToSave = [o.idsToSave,numel(o.allTimes)];            
        else
            %Use the number of steps specified (ignore the
            %timeStepSaveMode setting for restarts)
            nToAdd      = o.nSaveSteps - nOldSave;
            ids         = linspace(nOld,o.nTimeSteps,nToAdd+1);
            o.idsToSave = unique([o.idsToSave,round(ids(2:end))]);
        end
        o.nSaveSteps    = numel(o.idsToSave);
        o.times         = [o.times,zeros(1,o.nSaveSteps-nOldSave)];        
    else
        %New run -- Initialize variables for saving,
        %control and timing

        fOld            = oN.f(:,1);
        iSave           = 2;                
        firstStep       = 2;

        o.allTimes      = linspace(0,o.tMax,o.nTimeSteps);                                                    
        o.gmresFlags    = zeros(1,o.nTimeSteps);
        o.gmresRess     = zeros(1,o.nTimeSteps);
        o.gmresIters    = zeros(1,o.nTimeSteps);
        o.times         = zeros(1,o.nSaveSteps); 
    end
end
