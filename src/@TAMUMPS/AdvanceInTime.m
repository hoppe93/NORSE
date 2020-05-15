function AdvanceInTime(o,varargin)
    % Advance the system in time using the MUMPS direct solver.
    % Usage:
    %   AdvanceInTime()
    %   AdvanceInTime(true)
    %
    % If an argument is passed, it signifies a restart, which
    % affects array allocation and the time steps to perform.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    oN = o.norse;

    % Synch the settings with the NORSE object
    o.SynchSettings();

    % Initialize variables
    isARestart = (nargin == 2 && varargin{1});
    [fOld,iSave,firstStep] = o.InitializeTimeStepArrays(varargin{:});
    if isARestart
        tEndInv = oN.timing.matrixInversion;
        tEndAnalysis = oN.timing.mumpsAnalysis;
    else
        tEndInv = 0;
        tEndAnalysis = 0;
    end

    str = sprintf(' Beginning MUMPS time advance with %d steps\n   Step: ', o.nTimeSteps);
    o.PrintInfo(isARestart, str);

    %%% Advance in time %%%
    tStartLoop = tic;
    for iteration = firstStep:o.nTimeSteps
        t = o.allTimes(iteration);
        tOld = o.allTimes(iteration-1);
        isSaveStep = (iteration == o.idsToSave(iSave));
        o.iTimeStep = iSave;

        o.PrintProgress(isSaveStep, iSave, iteration);
        o.UpdateKineticEquation(fOld, t, tOld);

        tStartInv = tic;

        % MUMPS inversion
        % :: o.matrix * fOld = o.rhs

        % 1. Matrix analysis
        o.MUMPSAnalyze(o.matrix);
        tEndAnalysis = tEndAnalysis + toc(tStartInv);

        % 2. Factorization & 3. Solution
        if o.thresholdPivoting >= 0
            o.mumpsInfo.CNTL(1) = o.thresholdPivoting;
        end

        if o.precinv >= 0
            o.mumpsInfo.CNTL(2) = o.precinv;
        end

        fOld = o.MUMPSSolve(o.matrix, o.rhs);

        % %%%%%%%%%%%%%%%%%

        tEndInv   = tEndInv + toc(tStartInv);

        if isSaveStep  % Save the current state
            iSave = o.SaveState(iSave, fOld, t);
        end

        iSave = o.CheckAbortCriterion(fOld, t, isSaveStep, iSave, iteration);
        if o.abortFlag
            break
        end
    end

    o.dtsUsed = [o.dtsUsed, o.dt*ones(1, o.nTimeSteps-firstStep+1)];
    o.SynchAfterCalculation();

    % Save timings
    oN.timing.timeAdvance     = oN.timing.timeAdvance + toc(tStartLoop);
    oN.timing.matrixInversion = tEndInv;
    oN.timing.mumpsAnalysis   = tEndAnalysis;

    oN.Print('\n   Done!\n');

end
