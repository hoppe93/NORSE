classdef TAAdaptive < TimeAdvance
    % TAAdaptive -- Uses an adaptive time step based on the performance of
    %               GMRES to advance the NORSE system in time.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Usage: 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Initialization:
    %       TAAdaptive() -- empty object
    %
    %   In NORSE:
    %       TAAdaptive(oNORSE) 
    %       AdvanceInTime()
    %     
    %   
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Interface methods and output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods 
        AdvanceInTime(o,varargin)
            % Advance the system in time using indirect (GMRES) matrix
            % factorization with and adaptive time step. The time step is
            % dynamically adapted based on the number of iterations gmres
            % requires for convergence.
        
    end
    
    %%% Internal properties and methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        GMRESTolerance
        nStepsBetweenFactorizations
        dtIncreaseLimit         
        
        %LU factors
        L
        U
        P
        Q
        
        convergenceFlag = false
    end
    
    methods 
        function o = TAAdaptive(varargin)
            % Constructor.
            %
            % Usage:
            %   TAAdaptive()       -- empty object
            %   TAAdaptive(oNORSE) -- call to use in NORSE 
            
            o@TimeAdvance(varargin{:});
            o.CalculateNumberOfSteps();
        end
        
        CalculateNumberOfSteps(o)
            %Overloads the superclass method to account for that the number
            %of timesteps is unknown.
        SynchSettings(o)
            %Overloads the superclass method -- some additional properties
            %need to be specified
        [fOld,iSave,varargout] = InitializeTimeStepArrays(o,varargin)
            %Overloads the superclass method to account for that the number
            %of timesteps is unknown.
        dt = AdaptTimeStep(o,dt,nIterations)
            % Adapts time step based on the number of GMRES steps to
            % convergence. 
        ExtendSaveArrays(o,stepCounter,nStepGuess)
            % Allocates additional space for the save arrays. This is
            % needed in the adaptive time-step scheme, since it is not
            % known in advance how many time steps will be needed.
        TrimSaveArrays(o,stepCounter,iSave)
            %Overloads the superclass method
        f = PreconditionerFunc(o,rhs)
            % Efficiently computes the preconditioner used in the iterative
            % time advance schemes. Significantly speeds up the computation,
            % compared to supplying the L & U matrices directly.
        PrintConvergenceWarning(o,flag)
            % Prints a warning if the GMRES iteration raises an error flag.
        CleanUp(o)
            % Clear large internal variables that are not needed when
            % saving the NORSE object.
    end 
    
end
