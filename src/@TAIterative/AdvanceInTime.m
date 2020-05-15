function AdvanceInTime(o,varargin)
    % Advance the system in time using indirect (GMRES) matrix
    % factorization. The method converges in just a few iterations
    % Periodically, the system is solved exactly to provide a good
    % preconditioner. 
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
    isFactorizationStep = zeros(1,o.nTimeSteps);
    isFactorizationStep(firstStep:o.nStepsBetweenFactorizations:end) = 1;
    if isARestart
        %Restart, keep adding to previous timings
        tEndFac   = oN.timing.matrixFactorization; 
        tEndSolve = oN.timing.GMRES;                
    else
        %New run                
        tEndFac   = 0;
        tEndSolve = 0;                
    end            

    %Print info and prepare "live" plot            
    str = sprintf(' Beginning iterative time advance with %d steps\n   Step: ',o.nTimeSteps);
    o.PrintInfo(isARestart,str);

    %%% Advance in time %%%
    %Turn off a common warning, so that we can control when it is
    %shown using the GMRES output:
    warning off MATLAB:gmres:tooSmallTolerance
    maxit       = 25; %Maximum number of GMRES iterations
    tStartLoop  = tic; 
    for iteration = firstStep:o.nTimeSteps  
        t = o.allTimes(iteration);
        tOld = o.allTimes(iteration-1);
        isSaveStep = (iteration == o.idsToSave(iSave));
        o.iTimeStep = iSave;

        o.PrintProgress(isSaveStep,iSave,iteration);
        o.UpdateKineticEquation(fOld,t,tOld); %Build the system to solve

        %Periodically LU factorize the matrix, for use as a
        %preconditioner in later time steps, otherwise step in time 
        %using gmres
        if isFactorizationStep(iteration)
            tStartFac = tic;         
            [o.L,o.U,o.P,o.Q] = lu(o.matrix); 
            fOld = o.PreconditionerFunc(o.rhs);
            tEndFac = tEndFac+toc(tStartFac);
        else
            tStartSolve = tic;
            %Use L and U as preconditioners and the previous f as an initial guess
            [fOld,flag,o.gmresRess(iteration),iter] = ...
                    gmres(o.matrix,o.rhs,[],o.GMRESTolerance,maxit,...
                          @o.PreconditionerFunc,[],fOld); 
            o.gmresFlags(iteration) = flag;
            o.gmresIters(iteration) = iter(2);
            tEndSolve = tEndSolve+toc(tStartSolve);
            o.PrintConvergenceWarning(flag);
            o.PrintIterationWarning(iter);
        end

        if isSaveStep %Save the current state 
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
    o.dtsUsed = [o.dtsUsed,o.dt*ones(1,o.nTimeSteps-firstStep+1)]; %If this is a new run, the first entry will be o.dtsUsed=0
    o.SynchAfterCalculation();

    %Save timings
    oN.timing.timeAdvance         = oN.timing.timeAdvance + toc(tStartLoop);
    oN.timing.MatrixFactorization = tEndFac;
    oN.timing.GMRES 			  = tEndSolve;  

    oN.Print('\n   Done!\n');
    warning on MATLAB:gmres:tooSmallTolerance
end
