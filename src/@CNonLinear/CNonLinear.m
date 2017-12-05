classdef CNonLinear < CollisionOperator
    % CNONLINEAR -- Implements the relativistic nonlinear electron-
    %               electron collision operator of Braams and Karney, as
    %               well as a simple electron-ion operator.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Usage: 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Initialize:
    %   CNonLinear() -- empty object
    %
    %   In NORSE:
    %       CNonLinear(oNORSE)
    %       InitializeParts()
    %
    %   This collision operator must be rebuilt in each time step, via a 
    %   call to
    %       Assemble(t).
    %   The operator is then accessible through the class property C.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Interface methods and output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    methods 
        InitializeParts(o) 
            % Overloads the superclass method. Calculates the prefactors
            % that multiply the different potentials in the various terms
            % in the collision operator.
        Assemble(o,t) 
            % Overloads the superclass method. Assembles the precalculated
            % collision operator parts into one operator matrix.
    end
    
    %%% Internal properties and methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties         
        % Prefactors for terms of C
        d2fdp2
        dfdp
        d2fdxi2
        dfdxi
        d2fdpdxi
        f        
        ionD2fdxi2
        ionDfdxi
        
        preFactorEE 
        preFactorEI        
    end
    
    methods 
        function o = CNonLinear(varargin)
            % Constructor.
            %
            % Usage:
            %   CollisionOperator() -- empty object
            %   CollisionOperator(oNORSE) -- call to use in NORSE 
            
            o@CollisionOperator(varargin{:});            
        end     
    end
end
