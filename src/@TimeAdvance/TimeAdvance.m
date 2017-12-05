classdef (Abstract) TimeAdvance < matlab.mixin.Copyable
    % TimeAdvance -- Abstract class that handles the advancement of the 
    %                NORSE system in time. Subclasses need to implement the
    %                abstract methods defined here to be usable.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Usage: 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Initialization:
    %       TimeAdvance() -- empty object
    %
    %   In NORSE:
    %       TimeAdvance(oNORSE) 
    %       AdvanceInTime()
    %     
    % If oNORSE already contains a TimeAdvance object, its properties are
    % transferred to the new object. This is useful when changing 
    % timeAdvanceMode when continuing a calculation.
    %   
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Interface methods and output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        dt
        tMax
        
        nTimeSteps
        nSaveSteps
        
        times
        allTimes        
        dtsUsed   = 0
        
        gmresFlags        
        gmresIters
        gmresRess
        
        abortFlag = false
    end
    
    methods (Abstract) 
        AdvanceInTime(o,varargin)
            %This method performs the time advance and must be implemented
            %by a subclass in order for it to work in NORSE
    end
    
    %%% Internal properties and methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties 
        %Objects
        norse        
        
        %Time step arrays        
        idsToSave
        timesToSave
        
        %Elements of the kinetic equation
        matrix
        rhs
    end
    
    methods 
        function o = TimeAdvance(varargin)
            % Constructor.
            %
            % Usage:
            %   TimeAdvance() -- empty object
            %   TimeAdvance(oNORSE) -- call to use in NORSE            
            
            switch nargin
                case 0
                    %Nothing to do
                case 1
                    if isa(varargin{1},'NORSE')
                        o.norse = varargin{1};                        
                        switch o.norse.timeStepSaveMode
                            case {0,1}
                                %OK -- nothing to do
                            otherwise
                                error('Invalid time-step save mode.');
                        end
                        if ~isempty(o.norse.timeAdvance)
                            %Transfer information from an existing 
                            %timeAdvance object (needed to switch between
                            %time-advance modes at restarts)
                            o.norse.timeAdvance.CleanUp();
                            props = properties('TimeAdvance');
                            oTA = o.norse.timeAdvance;
                            for i = 1:numel(props)
                                o.(props{i}) = oTA.(props{i}); 
                            end                            
                        end
                    else
                        error('The first argument must be a NORSE object.'); 
                    end
                otherwise
                    error('Invalid number of input arguments');
            end
        
        end
        
        CalculateNumberOfSteps(o)
            
        SynchSettings(o)
            %Makes sure the settings are consistent with those in the
            %NORSE object.
        [fOld,iSave,varargout] = InitializeTimeStepArrays(o,varargin)
            %Initializes the time-step vector and various save arrays. Also
            %handles restarts.
        PrintInfo(o,isARestart,infoStr)
            % Prints the preliminary info about the time advance scheme.        
        PrintProgress(o,isSaveStep,iSave,iteration)
            % Prints the progress of the time advance at regular intervals.
        UpdateKineticEquation(o,fOld,t,tOld)
            % Updates and builds the matrix and right-hand side describing
            % the kinetic equation.
        iSave = SaveState(o,iSave,fOld,t)
            % Saves the state of the calculation and optionally plots the
            % parallel distribution.
        SaveStepData(o,iSave,f,t)
            % Saves the state of the solution at a given time step for
            % post-processing and visualization.
        iSave = CheckAbortCriterion(o,fOld,t,isSaveStep,...
                                                        iSave,iteration)
            % Check whether to abort the calculation because the system has
            % entered the slide-away regime.
        TrimSaveArrays(o,stepCounter,iSave)
            % Removes unused space at the end of the vectors/matrices used 
            % for saving time step data. This is used with the adaptive
            % time-step scheme or when aborting the calculation prematurely.
        SynchAfterCalculation(o)
            % Updates some fields in the NORSE object to reflect the final
            % state of the time advance.
        CleanUp(o)
            % Clear large internal variables that are not needed when
            % saving the NORSE object.       
        GetSize(o)
            % Calculates the size in memory of the object and prints it to
            % the console.
    end 
    
end
