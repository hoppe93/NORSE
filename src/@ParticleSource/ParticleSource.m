classdef ParticleSource < matlab.mixin.Copyable
    % PARTICLESOURCE -- Class that describes the particle source/sink in
    %                   NORSE needed to perform changes to the electron
    %                   density.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Usage: 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % Initialize with one of the calls:
    %   ParticleSource()       -- empty object
    %   ParticleSource(oNORSE) -- used in a NORSE calculation. Builds the
    %                             source with appropriate properties.
    %
    % Before starting a NORSE calculation, make sure Initialize() is
    % called to ensure consistency with NORSE parameters. 
    %
    % In each NORSE time step, call Update(t,tOld). The source operator
    % (with the appropriate magnitude) to include on the right-hand side of
    % the kinetic equation can then be accessed from the property source.
    %
    % The source magnitude is saved in the property densityChangeMagnitude, 
    % which is updated in every time step.        
    %    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Interface methods and output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties 
        source        
        densityChangeMagnitude = 0
    end
        
    methods 
        Initialize(o)
            % Initializes the quantities necessary for the particle
            % source/sink.
        Update(o,t,tOld)
            % Updates and returns the properly scaled sink.
    end    
    
    %%% Internal properties and methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties 
        %Objects
        norse
        grid
                            
        %%%%%%%%%%%%%%%%%
        
        particleSourceOperator        
        energyContentWBar
        particleSourceAP
        particleSourceDensityMoment
        
        %%%%%%%%%%%%%%%%%
        
        %Misc.
        gammaMin1
    end
    
    methods 
        function o = ParticleSource(varargin)
            % Constructor.
            %
            % Usage:
            %   ParticleSource()       -- empty object
            %   ParticleSource(oNORSE) -- call to use in NORSE, automatically
            %                             initializes the particle source
                        
            switch nargin
                case 0
                    %Nothing to do
                case 1
                    if isa(varargin{1},'NORSE')
                        o.norse = varargin{1};
                        o.grid  = o.norse.grid;                        
                    else
                        error('The first argument must be a NORSE object.'); 
                    end                    
                    o.Initialize();
                otherwise
                    error('Invalid number of input arguments');
            end            
        end
        
        BuildAndCalculateMagnitude(o,t,tOld)
            % Builds the particle source operator and calculates the
            % appropriate magnitude to perform specified changes to the
            % electron density.
        GetSize(o) 
            % Calculates the size in memory of the PARTICLESOURCE object by
            % looping through the fields and summing up the variable sizes,
            % and prints it to the console.
    end
end
