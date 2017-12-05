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
    
    %Initialize the arrays specified in the superclass
    [fOld,iSave,firstStep] = InitializeTimeStepArrays@TimeAdvance(o,varargin{:});
    
    %Also initialize some arrays specific to the inductive scheme
    o.inductiveCoefficients = struct('R1',[],'R2',[],'R3',[],...
                                         'R4',[],'R5',[]);                                         
    if nargin == 2 && varargin{1}
	%Restart, extend existing arrays
        oN.EHatInductive = [oN.EHatInductive,zeros(1,o.nSaveSteps-iSave+1)]; 
    else
	%New run -- initialize arrays from scratch
        oN.EHatInductive = zeros(1,o.nSaveSteps);
        oN.EHatInductive(1) = oN.EHat(o.allTimes(1));
    end
    if oN.LHat.isScalar
        LString = sprintf('LHat = %.1g',oN.LHat(0)); 
    else
        LString = 'a time-varying LHat';
    end
    oN.Print('   Calculating the electric field inductively with %s.\n',LString);
end
