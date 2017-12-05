classdef CollisionOperator < matlab.mixin.Copyable
    % COLLISIONOPERATOR -- Abstract class that implements the
    %                      electron-electron and electron-ion collision
    %                      operators in NORSE, as well as artificial
    %                      damping at the boundary at pMax.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Usage: 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Initialize:
    %   CollisionOperator() -- empty object
    %
    %   In NORSE:
    %       CollisionOperator(oNORSE)
    %       InitializeParts()
    %		SynchSettings()
    %
    %   The collision operator at time t is built via a call to
    %       Assemble(t).
    %   It is then accessible through the class property C.
    %
    %   The artificial dampening at pMax is influenced by the NORSE
    %   properties dampeningStrength and dampeningWidth. See NORSE.m for 
    %   details.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Interface methods and output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    properties 
        C
    end
    
    methods (Abstract)
        %These methods must be implemented in subclasses for them to be
        %usable in NORSE
        InitializeParts(o) 
            % Calculates the prefactors that multiply the different
            % potentials in the various terms in the collision operator. 
        Assemble(o,t) 
            % Assembles the precalculated collision operator parts into one
            % operator matrix. Uses the potentials calculated from the
            % distribution.        
    end
    
    methods
        function SynchSettings(o)
            % Makes sure settings are consistent with those in NORSE.
            % Should be called before the calculation starts (or at restart).
            
            %Nothing to do here in general
        end
    end
    
    %%% Internal properties and methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties 
        %Objects
        norse
        grid
        potentials
                      
        dampeningMat
    end
    
    methods 
        function o = CollisionOperator(varargin)
            % Constructor.
            %
            % Usage:
            %   CollisionOperator() -- empty object
            %   CollisionOperator(oNORSE) -- call to use in NORSE 
            
            switch nargin
                case 0
                    %Nothing to do
                case 1
                    if isa(varargin{1},'NORSE')
                        o.norse = varargin{1};
                        o.grid = o.norse.grid;
                        o.potentials = o.norse.potentials;
                        
                        o.InitializeArtificialDamping();
                    else
                        error('The first argument must be a NORSE object.'); 
                    end
                otherwise
                    error('Invalid number of input arguments');
            end
        end
        
        InitializeArtificialDamping(o)
            % Perpare an artifical dampening term for the boundary
            % condition at p_max. Can be used to reduce ringing in the tail
            % of the distribution.
        GetSize(o) 
            % Calculates the size in memory of the COLLISIONOPERATOR object
            % by looping through the fields and summing up the variable
            % sizes, and prints it to the console.            
    end
end
