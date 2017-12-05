classdef TALinearDirect < TimeAdvanceLinear
    % TALinearDirect -- Uses a direct matrix solve to advance the NORSE
    %                   system in time, using a linear collision operator.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Usage: 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Initialization:
    %       TALinearDirect() -- empty object
    %
    %   In NORSE:
    %       TALinearDirect(oNORSE) 
    %       AdvanceInTime()
    %     
    %   
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Interface methods and output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    methods 
        AdvanceInTime(o,varargin)        
            % Advance the system in time using direct (\) matrix
            % factorization.
    end
    
    %%% Internal properties and methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    methods 
        function o = TALinearDirect(varargin)
            % Constructor.
            %
            % Usage:
            %   TALinearDirect()       -- empty object
            %   TALinearDirect(oNORSE) -- call to use in NORSE 
            
            o@TimeAdvanceLinear(varargin{:});
            o.CalculateNumberOfSteps();
        end        
    end 
    
end
