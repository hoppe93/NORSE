classdef TimeAdvanceLinear < TimeAdvance
    % TimeAdvanceLinear -- Abstract class that handles the advancement of
    %                      the NORSE system in time when the linear
    %                      collision operator is used. Subclasses need to
    %                      implement the abstract method AdvanceInTime().
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Usage: 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Initialization:
    %       TimeAdvanceLinear() -- empty object
    %
    %   In NORSE:
    %       TimeAdvanceLinear(oNORSE) 
    %       TimeAdvanceLinear()
    %     
    %   
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Interface methods and output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%% Internal properties and methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    methods 
        function o = TimeAdvanceLinear(varargin)
            % Constructor.
            %
            % Usage:
            %   TimeAdvanceLinear()       -- empty object
            %   TimeAdvanceLinear(oNORSE) -- call to use in NORSE 
            
            o@TimeAdvance(varargin{:});            
        end
        
        UpdateMatrix(o,fOld,t,tOld)
            %Rebuilds the matrix describing the E-field, synchrotron and
            %collision operator terms. Only needs to be called when
            %parameters change.
        UpdateRHS(o,fOld,t,tOld)
            %Updates the right-hand side of the kinetic equation, which
            %changes in every time step.
        parsHaveChanged = FindWhenParametersHaveChanged(o,firstStep)
            %Evaluates when the time-dependent physical parameters have
            %changed so that the matrix will need to be rebuilt.
    end 
    
end
