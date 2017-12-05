function AdvanceInTime(o,varargin)
    % Advance the system describing F and EHat in time using Newton's
    % method. This method calculates the electric field self-consistently
    % together with the distribution. 
    %
    % Usage:
    %   AdvanceInTime()
    %   AdvanceInTime(true)
    %
    % If an argument is passed, it signifies a restart, which
    % affects array allocation and the time steps to perform.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    oN = o.norse;
    
    %Synch the settings with the NORSE object
    o.SynchSettings();

    %Initialize variables
    isARestart = (nargin == 2 && varargin{1});
    [fOld,iSave,firstStep] = o.InitializeTimeStepArrays(varargin{:});
    if isARestart                
        tEndInv = oN.timing.matrixInversion;
    else                
        tEndInv = 0;            
    end

    %Print info and prepare "live" plot            
    str = sprintf(' Beginning direct time advance with %d steps\n   Step: ',o.nTimeSteps);
    o.PrintInfo(isARestart,str);

    %%% Advance in time %%%
    EHatInd = oN.EHatInductive(iSave-1);
    tStartLoop = tic;    
    for iteration = firstStep:o.nTimeSteps    
        t = o.allTimes(iteration);
        tOld = o.allTimes(iteration-1);
        isSaveStep = (iteration == o.idsToSave(iSave));

        o.PrintProgress(isSaveStep,iSave,iteration); 
        o.UpdateKineticEquation(fOld,t,tOld); 
                            %Build the components of the system to solve
        
        %Take a time step by iterating using Newton's method        
        tStartInv = tic;        
        for iNewton = 1:o.nNewtonSteps
            [fOld,EHatInd] = o.TakeNewtonStep(fOld,EHatInd);
        end
        tEndInv   = tEndInv+toc(tStartInv);

        if isSaveStep %Save the current state 
            oN.EHatInductive(iSave) = EHatInd;
            iSave = o.SaveState(iSave,fOld,t);            
        end

        %Check whether to abort the calculation
        iSave = o.CheckAbortCriterion(fOld,t,isSaveStep,iSave,iteration);
        if o.abortFlag 
            break
        end 
    end

    %Save the time step lengths used. This is for compliance with
    %the adaptive time-advance scheme.
    o.dtsUsed = [o.dtsUsed,o.dt*ones(1,o.nTimeSteps-firstStep+1)];
    o.SynchAfterCalculation();

    %Save timings
    oN.timing.timeAdvance     = oN.timing.timeAdvance + toc(tStartLoop);
    oN.timing.matrixInversion = tEndInv;                        

    oN.Print('\n   Done!\n');
end
