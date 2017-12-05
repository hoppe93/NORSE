classdef Potentials < matlab.mixin.Copyable
    % POTENTIALS -- Handles the calculation of the potentials of the
    %               Braams & Karney collision operator in NORSE.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Usage: 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Initialization:
    %       Potentials() -- empty object
    %
    %   In NORSE:
    %       Potentials(oNORSE) 
    %       GeneratePotentialMatrices()
    %     
    % The potentials must be updated in each time step via a call to 
    %   Update(fls).
    % They are then accessible through the properties
    %   Upsilon0, Upsilon1, Upsilon2, UpsilonPlus, UpsilonMinus & Pi.
    %   
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Potentials and interface methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties 
        Upsilon0
        Upsilon1
        Upsilon2
        
        UpsilonPlus
        UpsilonMinus
        Pi
    end
    
    methods 
        GeneratePotentialMatrices(o)
            % Generates the matrices needed to efficiently calculate the
            % potentials from the distribution.
        Update(o,fls) 
            % Updates the potentials based on the distribution fls in the
            % finite-difference--Legendre-mode basis.
    end
    
    %%% Internal properties and methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties 
        %Objects
        norse
        grid
                
        %Matrices for calculating the potentials from the distribution
        Upsilon0Matrix
        Upsilon1Matrix
        Upsilon2Matrix
        PiMatrix
        
        %Misc. for building calculation matrices
        jAndYMatrices
        La_lBase
    end
    
    methods 
        function o = Potentials(varargin)
            % Constructor.
            %
            % Usage:
            %   Potentials() -- empty object
            %   Potentials(oNORSE) -- call to use in NORSE 
            
            switch nargin
                case 0
                    %Nothing to do
                case 1
                    if isa(varargin{1},'NORSE')
                        o.norse = varargin{1};
                        o.grid = o.norse.grid;
                    else
                        error('The first argument must be a NORSE object.'); 
                    end
                otherwise
                    error('Invalid number of input arguments');
            end
        end
        
        sPot = GeneratePotentialMatricesForL(o,l)
            % Calculates the matrices needed to calculate the various
            % potentials from the distribution for a given Legendre mode l.            
        s = GenerateBoundaryConditionsForInt(o,l)
            % Uses Eq. (31) in Braams & Karney [PoF B, 1 (1989), 1355] to
            % calculate boundary conditions for the differential
            % calculation of the potentials at p=pMax. 
        NLA = GetNLABoundary(o,l,a) 
            % Generate a vector NLA with the p' dependence of N_l,a at pMax.
            % This can then be used to integrate over p'.
        CalculateJAndY(o)
            % Calculate arrays containing j_l[k]* and y_l[k]* for all grid
            % points and all Legendre modes.        
        m = PrepareBoundaryConds(o,m,l)
            % Sets the proper l-dependent boundary conditions at p=0 and
            % p=pmax in the matrix m. l is the Legendre-mode index.
        GetSize(o) 
            % Calculates the size in memory of the POTENTIALS object by
            % looping through the fields and summing up the variable sizes,
            % and prints it to the console.
    end 
    
end
