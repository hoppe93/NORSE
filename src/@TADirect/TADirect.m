classdef TADirect < TimeAdvance
    % TADirect -- Uses a direct matrix solve to advance the NORSE system in
    %             time.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Usage: 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Initialization:
    %       TADirect() -- empty object
    %
    %   In NORSE:
    %       TADirect(oNORSE) 
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

        function o = TADirect(varargin)
            % Constructor.
            %
            % Usage:
            %   TADirect()       -- empty object
            %   TADirect(oNORSE) -- call to use in NORSE 
            
            o@TimeAdvance(varargin{:});
            o.CalculateNumberOfSteps();
        end        

        InitializeMUMPS(o, mumpsPath)

    end 
    
end
