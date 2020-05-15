classdef HeatSink < matlab.mixin.Copyable
    % HEATSINK -- Class that describes the heat sink in NORSE and related
    %             quantities such as the contribution to the source
    %             magnitude from various effects.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Usage: 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % Initialize with one of the calls:
    %   HeatSink()       -- empty object
    %   HeatSink(oNORSE) -- used in a NORSE calculation. Builds the heat  
    %                       sink with appropriate properties.
    %
    % Before starting a NORSE calculation, make sure SynchSettings() is
    % called to ensure consistency with NORSE parameters. 
    %
    % In each NORSE time step, call Update(f,t,tOld). The sink operator
    % (with the appropriate magnitude) to include in the kinetic equation
    % can then be accessed from the property sink.
    %
    % The struct describing the sink magnitude, as well as the
    % contribution to it from different effects, is saved in the
    % property energyChangeMagnitude, which is updated in every time step.
    %
    % See NORSE.m for a description of the settings affecting the heat
    % sink.
    %    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Interface methods and output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties 
        sink        
        energyChangeMagnitude
    end
        
    methods
        SynchSettings(o)
            % Makes sure settings are consistent with those in NORSE.
            % Should be called before the calculation starts (or at restart).
        Update(o,f,t,tOld)
            % Updates and returns the properly scaled sink.    
    end    
    
    %%% Internal properties and methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties 
        %Objects
        norse
        grid
        particleSource        
                            
        %%%%%%%%%%%%%%%%%
        
        %Flags
        includeHeatSink
        includeTempChanges
        calculateDensityRelatedChanges
        includeDensityRelatedChanges
        enforceStrictHeatConservation
        restrictHeatSinkRate
        
        %%%%%%%%%%%%%%%%%
        
        heatSinkOperator
        preFactor
        energyChangeInt
        densityInt
        magnitudeStructTemplate
        rateNormalization
        rateCutOff
        
    end
    
    methods 
        function o = HeatSink(varargin)
            % Constructor.
            %
            % Usage:
            %   HeatSink()       -- empty object
            %   HeatSink(oNORSE) -- call to use in NORSE, automatically
            %                       builds the sink operator
                        
            switch nargin
                case 0
                    %Nothing to do
                case 1
                    if isa(varargin{1},'NORSE')
                        o.norse          = varargin{1};
                        o.grid           = o.norse.grid;
                        o.particleSource = o.norse.particleSource;                        
                    else
                        error('The first argument must be a NORSE object.'); 
                    end
                    o.SynchSettings();
                    o.Assemble();
                otherwise
                    error('Invalid number of input arguments');
            end
            %Create a template struct for saving the source magnitude and
            %its components
            o.magnitudeStructTemplate = struct('sink',0,'E',0,'C',0,'synch',0,...
                                          'tempChange',0,'densityChange',0,...
                                          'correction',0,'restriction',0,...
                                          'normalization',0);
        end
        
        Assemble(o,varargin)
            % Builds the operator describing the heat sink, as well as
            % quantities needed to determine the heat sink magnitude based 
            % on the distribution and forces or parameter changes.                    
        CalculateEnergyChangeMagnitude(o,f,t,tOld)
            % Calculates the appropriate magnitude of the heat sink to
            % counteract the Ohmic heating, synchrotron losses, and effects
            % of collisions, as well as carry out temperature or
            % density-change-related energy changes.            
        GetSize(o) 
            % Calculates the size in memory of the HEATSINK object by looping
            % through the fields and summing up the variable sizes, and
            % prints it to the console.         
    end
end
