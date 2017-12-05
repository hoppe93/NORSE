function PrintProgress(o,isSaveStep,iSave,iteration)
    % Prints the progress of the time advance at regular intervals.
    %
    % Usage:
    %   PrintProgress(isSaveStep,iSave,iteration)
    %
    % isSaveStep is a boolean determining whether to print
    % something, iSave is an index used to introduce a line-break
    % every 10 printouts, and iteration is the iteration number.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isSaveStep 
        o.norse.Print('%d...',iteration);        
        if mod(iSave-1,10) == 0
            o.norse.Print('\n         ')
        end
    end
end
