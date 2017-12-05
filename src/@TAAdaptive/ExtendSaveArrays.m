function ExtendSaveArrays(o,stepCounter,nStepGuess)
    % Allocates additional space for the save arrays. This is
    % needed in the adaptive time-step scheme, since it is not
    % known in advance how many time steps will be needed.
    %
    % Usage:
    %   ExtendSaveArrays(stepCounter,nStepGuess)
    %
    % stepCounter is the total number of steps taken so far and
    % nStepGuess is the initial guess for the number of steps
    % needed.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Every time we have filled up the empty entries, add an
    %additional nStepGuess entries
    if ~mod(stepCounter,nStepGuess)                 
        o.gmresFlags  = [o.gmresFlags,zeros(1,nStepGuess)];
        o.gmresRess   = [o.gmresRess,zeros(1,nStepGuess)];
        o.gmresIters  = [o.gmresIters,zeros(1,nStepGuess)];
        o.dtsUsed     = [o.dtsUsed,zeros(1,nStepGuess)];
        o.allTimes    = [o.allTimes,zeros(1,nStepGuess)];
    end 
end
