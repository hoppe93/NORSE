function [fOld,iSave,tVars] = InitializeTimeStepArrays(o,varargin)
    %Initializes the time-step vector and various save arrays for
    %the adaptive time-advancement scheme. Also handles restarts.
    %
    % Usage: 
    %   [fOld,iSave,tVars] = InitializeTimeStepArrays()
    %           -- The struct tVars contains several fields that
    %              are related to time, such as counters, and dt
    %                [...] = o.InitializeTimeStepArrays(true) -- restart
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    oN = o.norse;
    tVars = struct();
    tVars.nStepGuess = 100;
    if nargin == 2 && varargin{1}
        %Restart -- extend arrays rather than initialize them                        
        if isempty(o.timesToSave)
            %A different time-advancement scheme was used
            %before the restart. Initialize an array of the
            %saved times for consistency.
            o.timesToSave = o.times;
        end

        %Extract information from the existing save arrays
        fOld    = oN.f(:,end);
        iSave   = numel(o.times)+1;
        nOld    = numel(o.allTimes);
        tOld    = o.times(end);

        %Determine how many and what steps to save
        if o.nSaveSteps < nOld
            %Save only the final time (and keep all the old
            %saved steps)
            o.timesToSave = [o.timesToSave,o.tMax];
            nToAdd = 1;
        else
            %Use the number of steps specified (ignore the
            %timeStepSaveMode setting for restarts)
            nToAdd = o.nSaveSteps - nOld;                            
            timesToAdd    = linspace(tOld,o.tMax,nToAdd+1);
            o.timesToSave = [o.timesToSave,timesToAdd(2:end)];                            
        end                        

        %Set the various times and counters
        if o.dtsUsed(2) ~= o.dt
            %The dt has been changed by the user. Use the
            %supplied value
            tVars.dt      = o.dt;
        else
            %Start with the dt used at the end, not the 
            %initial dt
            tVars.dt      = o.dtsUsed(end); 
        end
        tVars.stepCounter = nOld+1;                         
        tVars.tOld        = tOld;                        
        tVars.t           = tOld+tVars.dt;
        tVars.iter        = [0,o.gmresIters(end)];
        tVars.factorizationCounter = o.nStepsBetweenFactorizations;%Triggers a factorization in the first iteration                                   
        tVars.tEndFac     = oN.timing.matrixFactorization;
        tVars.tEndSolve   = oN.timing.GMRES;

        %Extend existing arrays (guess the number of time
        %steps needed).
        nG = tVars.nStepGuess;
        o.gmresFlags  = [o.gmresFlags,zeros(1,nG)];
        o.gmresRess   = [o.gmresRess,zeros(1,nG)];
        o.gmresIters  = [o.gmresIters,zeros(1,nG)];
        o.dtsUsed     = [o.dtsUsed,zeros(1,nG)];
        o.allTimes    = [o.allTimes,zeros(1,nG)];
        o.times       = [o.times,zeros(1,nToAdd)]; 
    else                
        %Initialize variables for saving, control and timing                        
        fOld          = oN.f(:,1);            
        iSave         = 2;

        %Set the various times and counters
        tVars.factorizationCounter = 0;                            
        tVars.stepCounter   = 2;              
        tVars.iter          = [0,0];
        tVars.dt            = o.dt;
        tVars.t             = o.dt;
        tVars.tOld          = 0;
        tVars.tEndFac       = 0;
        tVars.tEndSolve     = 0;

        %Initialize arrays (guess the number of time steps
        %needed).
        nG = tVars.nStepGuess;
        o.gmresFlags  = zeros(1,nG);
        o.gmresRess   = zeros(1,nG);
        o.gmresIters  = zeros(1,nG);
        o.dtsUsed     = zeros(1,nG);
        o.allTimes    = zeros(1,nG);

        switch oN.timeStepSaveMode
            case 0 %Linear
                o.timesToSave = linspace(0,o.tMax,o.nSaveSteps);
            case 1 %Logarithmically increasing
                o.timesToSave = logspace(log10(o.dt),log10(o.tMax),o.nSaveSteps);
                o.timesToSave(1) = 0;
            otherwise
                error('Invalid timeStepSaveMode.');
        end                         
    end       
end
