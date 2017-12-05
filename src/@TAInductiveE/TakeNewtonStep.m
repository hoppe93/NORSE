function [fOld,EHatInd] = TakeNewtonStep(o,fOld,EHatInd)
    % Builds and solves the necessary matrices to perform one Newton
    % iteration for the determination of the coupled system of f and the
    % inductive electric field.
    %
    % Usage: 
    %   [fOld,EHatInd] = TakeNewtonStep(o,fOld,EHatInd)
    %   
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


    mS = o.norse.grid.matSize;
    
    %Build the residual vector
    R = [ (o.inductiveCoefficients.R1 ...
           + EHatInd*o.inductiveCoefficients.R2)*fOld ...
         + o.inductiveCoefficients.R3;...
         ...
         EHatInd + o.inductiveCoefficients.R4*fOld ...
         + o.inductiveCoefficients.R5];
     
    %Build the Jacobian matrix 
    D = sparse(mS+1,mS+1);
    D(1:mS,1:mS) = o.inductiveCoefficients.R1 ...
                 + EHatInd*o.inductiveCoefficients.R2;
    D(1:mS,end)  = o.inductiveCoefficients.R2*fOld;
    D(end,1:mS)  = o.inductiveCoefficients.R4;
    D(end,end)   = 1;
    
    %Solve for the updated state
    state   = [fOld;EHatInd] - D\R;
    fOld    = state(1:end-1);
    EHatInd = state(end);    
end
