classdef CLinear < CollisionOperator
    % CLINEAR -- Implements a linearized version of the Braams and Karney
    %            electron-electron collision operators (due to Pike and
    %            Rose) as well as a simple electron-ion operator. The e-e
    %            operator is relativistic and conservative, and does not
    %            need to be updated in each time step (unless the
    %            temperature or density changes).
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Usage: 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Initialize:
    %   CLinear() -- empty object
    %
    %   In NORSE:
    %       CLinear(oNORSE)
    %       InitializeParts()
    %		SynchSettings()
    %
    %   To build the collision operator, call
    %       Assemble(t).
    %   The operator is then accessible through the class property C.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Interface methods and output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    properties %Settings -- see description in the NORSE object comments
        potentialCutOff           = 1e-30
        specialFunctionEvaluation = 1
        includeFPTerm             = 1
        includeTPTerm             = 1
    end
        
    methods 
        InitializeParts(o) 
            % Overloads the superclass method. Calculates the prefactors
            % that multiply the different potentials in the various terms
            % in the collision operator.
        Assemble(o,t) 
            % Overloads the superclass method. Assembles the precalculated
            % collision operator parts into one operator matrix.
        SynchSettings(o)
            % Makes sure the settings are consistent with those of the
            % NORSE object.
    end
    
    %%% Internal properties and methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties 
        % Grid constants for FP term of e-e operator
        dUpsilon0
        dUpsilon2
        dUpsilonMin
        dUpsilonPlus
        dPi
        
        % Grid constants for TP term of e-e operator
        d2fdxi2
        dfdxi
        
        % Grid constants for e-i operator
        ionD2fdxi2
        ionDfdxi
                
        %Prefactors
        preFactorTP
        preFactorFP 
        preFactorEI        
        
        % Matrices for calculating combined potentials
        UpsilonMinMatrix
        UpsilonPlusMatrix
        
        % Matrices for calculating weighted potentials 
        weightedU0Mat
        weightedU2Mat
        weightedUMinMat
        weightedUPlusMat
        weightedPiMat
        
        % Special functions for TP term
        GPara
        GPerp
        GK
        dGParaDP
        dGKDP
    end
    
    methods 
        function o = CLinear(varargin)
            % Constructor.
            %
            % Usage:
            %   CollisionOperator() -- empty object
            %   CollisionOperator(oNORSE) -- call to use in NORSE 
            
            o@CollisionOperator(varargin{:});            
        end
        
        InitializeTestParticleTerm(o)
            % Calculates the time-independent prefactors for the
            % test-particle term of the electron-electron operator.
        InitializeFieldParticleTerm(o)
            % Calculates the time-independent prefactors for the
            % field-particle term of the electron-electron operator.        
        CalculateSpecialFunctions(o,t)
            % Determines the special functions G_K, G_Para, G_Perp (as well
            % as the derivatives dG_K/dp and dG_Para/dp) needed for the
            % test-particle term of the linear e-e collision operator.
        CalculateWeightedPotentialMatrices(o,t)
            % Weighs the matrices used to calculate the potentials from f
            % in the Legendre basis by an exponential factor. This reduces
            % the number of nonzero elements in the field-particle term.
    end
    
    methods (Static)
        cumInt = CumSimpsonsRule(f,p)
        % Calculate the cumulative integral of a function f(p) over the
        % grid p.
    end
end
