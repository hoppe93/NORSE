classdef TAIterative < TimeAdvance
    % TAIterative -- Uses an iterative matrix solve to advance the NORSE
    %                system in time.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Usage: 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Initialization:
    %       TAIterative() -- empty object
    %
    %   In NORSE:
    %       TAIterative(oNORSE) 
    %       AdvanceInTime()
    %     
    %   
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Interface methods and output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods 
        AdvanceInTime(o,varargin)
            % Advance the system in time using indirect (GMRES) matrix
            % factorization. 
    end
    
    %%% Internal properties and methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties 
        GMRESTolerance
        nStepsBetweenFactorizations
        
        %LU factors
        L
        U
        P
        Q
        
        convergenceFlag = false
        iterationsFlag = false
    end
    
    methods 
        function o = TAIterative(varargin)
            % Constructor.
            %
            % Usage:
            %   TAIterative()       -- empty object
            %   TAIterative(oNORSE) -- call to use in NORSE 
            
            o@TimeAdvance(varargin{:});
            o.CalculateNumberOfSteps();
        end
        
        SynchSettings(o)
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
