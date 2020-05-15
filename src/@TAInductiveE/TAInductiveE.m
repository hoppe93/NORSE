classdef TAInductiveE < TimeAdvance
    % TAInductiveE -- Uses Newton's method to solve for f and the electric
    %                 field consistently.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Usage: 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Initialization:
    %       TAInductiveE() -- empty object
    %
    %   In NORSE:
    %       TAInductiveE(oNORSE) 
    %       AdvanceInTime()
    %     
    %   
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Interface methods and output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    methods 
        AdvanceInTime(o,varargin)
            % Advance the system in time by calculating the electric field
            % self-consistently together with the distribution by iterating
            % using Newton's method.        
    end
    
    %%% Internal properties and methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties         
        inductiveCoefficients
        nNewtonSteps
        
        % MUMPS settings
        mumpsInfo
        useMumps = 0
        precinv = -1            % Precision of inversion: < 0 means sqrt(eps)
        thresholdPivoting = -1  % Threshold for pivoting: < 0 means MUMPS decides
    end
    
    methods 
        function o = TAInductiveE(varargin)
            % Constructor.
            %
            % Usage:
            %   TAInductiveE()       -- empty object
            %   TAInductiveE(oNORSE) -- call to use in NORSE 
            
            o@TimeAdvance(varargin{:});
            o.CalculateNumberOfSteps();

            % Initialize MUMPS
            oN = varargin{1};
            o.InitializeMUMPS(oN.mumpsPath);
        end
        
        SynchSettings(o)
            %Overloads the superclass method
        [fOld,iSave,firstStep] = InitializeTimeStepArrays(o,varargin)
            %Overloads the superclass method             
        UpdateKineticEquation(o,fOld,t,tOld)
            % Updates and builds the matrix and right-hand side describing
            % the kinetic equation in the case including an inductive E.
        [fOld,EHatInd] = TakeNewtonStep(o,fOld,EHatInd)
            % Builds and solves the necessary matrices to perform one
            % Newton iteration for the determination of f and E.
            
        
        InitializeMUMPS(o,mumpsPath);
        % Initialize the MUMPS solver
        
        MUMPSAnalyze(o,matrix);
        % Run the MUMPS analysis stage on the matrix (required before
        % solve)
        
        f = MUMPSSolve(o,matrix,rhs);
        % Solve the linear system using MUMPS
    end 
    
end
