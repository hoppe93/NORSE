function AdvanceInTime(o,varargin)
    % Advance the system in time using indirect (GMRES) matrix
    % factorization with and adaptive time step. Periodically, the
    % system is solved exactly to provide a good preconditioner.
    % The time step is dynamically adapted based on the number of
    % iterations gmres requires for convergence.
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

    %Print info and prepare "live" plot
    isARestart = (nargin == 2 && varargin{1});
    str = sprintf(' Beginning iterative time advance with adaptive steps. \n   Step: ');
    o.PrintInfo(isARestart,str);            

    [fOld,iSave,sT] = o.InitializeTimeStepArrays(varargin{:});

    %%% Advance in time %%%
    % Here we do not know in advance the number of time steps to
    % take, which modifies the algorithm compared to the direct and
    % iterative schemes. The important quantity here is the time t,
    % rather than the step index.

    %Turn off a common warning, so that we can control when it is
    %shown using the GMRES output:
    warning off MATLAB:gmres:tooSmallTolerance
    maxit       = 25; %Maximum number of GMRES iterations            
    tStartLoop  = tic;                
    while sT.t <= (o.tMax + sT.dt)
        isSaveStep = (sT.t >= o.timesToSave(iSave));
        o.iTimeStep = iSave;

        o.PrintProgress(isSaveStep,iSave,sT.stepCounter);
        o.UpdateKineticEquation(fOld,sT.t,sT.tOld); %Build the system to solve
        o.dtsUsed(sT.stepCounter) = sT.dt; %Save the time step used

        %Periodically LU factorize the matrix, for use as a
        %preconditioner in later time steps, otherwise step in time 
        %using gmres                
        if (sT.factorizationCounter == o.nStepsBetweenFactorizations) ...
                                            || (sT.stepCounter == 2)
            tStartFac = tic;         
            [o.L,o.U,o.P,o.Q] = lu(o.matrix); 
            fOld = o.PreconditionerFunc(o.rhs);
            sT.tEndFac = sT.tEndFac+toc(tStartFac);
            sT.factorizationCounter = 0;                                           
            sT.dt = o.AdaptTimeStep(sT.dt,sT.iter(2)); 
                                %Adapt time step based on the
                                %number of iterations in the
                                %previous call to gmres
        else
            tStartSolve = tic;
            %Use L and U as preconditioners and the previous f as an initial guess
            [fOld,flag,o.gmresRess(sT.stepCounter),sT.iter] = ...
                    gmres(o.matrix,o.rhs,[],o.GMRESTolerance,maxit,...
                          @o.PreconditionerFunc,[],fOld); 
            o.gmresFlags(sT.stepCounter) = flag;
            o.gmresIters(sT.stepCounter) = sT.iter(2);
            sT.tEndSolve = sT.tEndSolve+toc(tStartSolve);                                        
            o.PrintConvergenceWarning(flag);
        end                

        if isSaveStep %Save the current state 
            o.idsToSave(iSave) = sT.stepCounter; 
                        %This is redundant, but required for
                        %cross-compatibility with non-adaptive
                        %time-advance schemes
            iSave = o.SaveState(iSave,fOld,sT.t);                    
        end

        %Prepare for next time step
        o.allTimes(sT.stepCounter) = sT.t;
        sT.stepCounter = sT.stepCounter+1;
        sT.factorizationCounter = sT.factorizationCounter+1;
        sT.tOld = sT.t;
        sT.t = sT.t+sT.dt;           
        o.ExtendSaveArrays(sT.stepCounter,sT.nStepGuess)   

        %Check whether to abort the calculation
        iSave = o.CheckAbortCriterion(fOld,sT.t,isSaveStep,...
                                                iSave,sT.stepCounter);
        if o.abortFlag
            break
        end                                
    end                        

    %Save timings
    oN.timing.timeAdvance 		  = oN.timing.timeAdvance + toc(tStartLoop);
    oN.timing.matrixFactorization = sT.tEndFac;
    oN.timing.GMRES 		 	  = sT.tEndSolve;

    o.TrimSaveArrays(sT.stepCounter,iSave); 
    o.SynchAfterCalculation();

    oN.Print('\n   Done!\n');
    warning on MATLAB:gmres:tooSmallTolerance
end      
