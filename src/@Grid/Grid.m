classdef Grid < matlab.mixin.Copyable
    % GRID -- Class that describes the finite-difference grid in a NORSE
    %         calculation, and associated quantities such as
    %         finite-difference derivatives. Also contains various methods
    %         for manipulating the grid.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Usage: 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Initialize an empty object:
    %   Grid() -- Specify settings and then call InitializeGrid() to
    %             construct the grid.
    % Call to use in NORSE:
    %   Grid(oNORSE)
    %
    % Initialize a grid with given parameters and standard settings:
    %   Grid(nP,nXi,pMax,pGridMode,xiGridMode) 
    %   Grid(nP,nXi,pMax,pGridMode,xiGridMode,pGridParameter)
    %
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Settings: 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % nP             -- Number of points in the p grid.
    % pMax           -- Maximum value of the p grid (p = gamma*v/c).
    % pGridMode      -- Setting for how the nP points between 0 and pMax 
    %                   are chosen. Default: 1
    %                    0: Uniform
    %                    1: Approximately quadratic increase with p 
    %                       (p = s^2 + pGridParameter*s, with s uniform)
    %                    2: Approximately cubic increase with p 
    %                       (p = s^3 + pGridParameter*s, with s uniform)
    %                    3: Approximately quartic increase with p 
    %                       (p = s^4 + pGridParameter*s, with s uniform)
    %                    4: A grid with a tanh step in spacing, giving a 
    %                       dense grid at low p and a sparse grid at high p
    % pGridParameter -- Parameter used in pGridMode = 1-3. Default: 0.2
    % pStepParams    -- Struct of parameters for pGridMode = 4, with fields:
    %   - stepLocation   -- The value of y at which the step in spacing is
    %                       located. y=gamma v/v_th is the thermal momentum.
    %                       Default: 5   
    %   - bulkSpacing    -- Determines the spacing in the bulk (in units of 
    %                       y). Default: 0.01
    %   - Theta          -- The temperature at which the above parameters
    %                       are defined. Default: 0.01 (5.11 kev)
    %   - stepSharpness  -- Determines the rate of change of the spacing at 
    %                       the step. Default: 0.1 (arb. units)
    %
    % nXi            -- Number of points in the xi grid (xi = cos(theta)).
    % xiGridMode     -- Specifies how the nXi grid points between -1 and 1 are
    %                   chosen. Default: 1
    %                    0: Uniform (this mapping requires that nXi is odd)
    %                    1: A sigmoid, denser close to xi=1 
    %                    2: Uniform in the polar angle (arccos(\xi))  
    %
    % stencil        -- Stencil to use for finite-difference derivatives:
    %                   3,5 or 7. Default: 5
    % quadrature     -- Quadrature to use for integration on the grid(s):
    %                   'simpson' or 'trapezoid'. Default: 'simpson'
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Grid properties: 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % A quantity f on the grid can be represented in two ways; on a 2D grid
    % of size [nP,nXi], so that
    %   f = [f(p,xi=-1),...,f(p,xi=1)],
    % (where p denotes all p between p_0=0 and p_nP=pMax), or in vector
    % representation where 
    %   f = [f(pp,xi_1);f(pp,xi_2);...;f(pp,xi_nXi);f(p0)],
    % and pp denotes all p between p_1 and p_nP=pMax (note the special
    % treatment of the point at the origin of the coordinate system). The
    % vector has dimensions [matSize,1].
    % 
    % Below, array sizes are shown in brackets.
    %
    % matSize  -- Number of elements in the grid ((nP-1)*nXi+1). The linear
    %             system solved in NORSE has size [matSize,matSize].
    %
    % %%% p grid %%%
    % p        -- Points in the p grid. [nP,1]
	% gamma    -- Relativistic mass factor (gamma = sqrt(p^2+1)). [nP,1]
    % dp       -- p(2)-p(1). [1,1]    
    % pWeights -- Quadrature weights used for performing integrals over p.
    %             [nP,1]
    % ddp      -- Finite-difference derivative on p grid. [nP,nP]
    % d2dp2    -- Finite-difference second derivative on p grid. [nP,nP]
    % p2D      -- Values of p on the entire 2D grid. [nP,nXi]
    % pBig     -- p2D mapped to a vector. [matSize,1]
    % gammaBig -- gamma on the entire grid. [matSize,1]    
    % idsPMax  -- Indices for the location of the points p=pMax for all xi 
    %             in pBig and other vectors on the same form. [1,nXi]
    %     
    % %%%%%%%%%%%%%%%%%
    %   
    % %%% xi grid %%%
    % xi        -- Points in the xi grid. [nXi,1]
    % xiWeights -- Quadrature weights used for performing integrals over xi.
    %              [nXi,1]
    % ddxi      -- Finite-difference derivative on xi grid. [nXi,nXi]
    % d2dxi2    -- Finite-difference second derivative on xi grid. [nXi,nXi]
    % xi2D      -- Values of xi on the entire 2D grid. [nP,nXi]
    % xiBig     -- xi2D mapped to a vector. [matSize,1]
    % xi0Id     -- Index for the location of xi=0 in xi. [1,1]
    %    
    % %%%%%%%%%%%%%%%%%%
    %           
    % %%% p_para, p_perp space %%%
    % pPara2D   -- The value of p_parallel = p*xi on the 2D grid. [nP,nXi]
    % pPerp2D   -- The value of p_perpendicular = p*sqrt(1-xi^2) on the 2D 
    %              grid. [nP,nXi]
    %
    % %%%%%%%%%%%%%%%%%%
    %    
    % %%% Derivative and integral operators %%%
    % 
    % The following matrices represent finite-difference derivatives on 
    % the entire grid (in vector representation). [matSize,matSize]
    % ddpMat     
    % ddxiMat    
    % d2dp2Mat   
    % d2dxi2Mat  
    % d2dpdxiMat 
    %
    % The following vectors represent quadrature weights for performing 
    % integrals over a vector f on the 2D grid. [matsize,1]
    %
    % intdp      -- 4*pi*int_0^pMax p^2 f dp    
    % intdpdxi   -- 2*pi*int_-1^1 int_0^pMax p^2 f dp dxi 
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Useful methods: 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % InitializeGrid() -- Builds the grid and all the associated quantities
    %                     above based on the settings specified in the
    %                     object properties.
    %
    % Visualize()      -- Plots the grid points and their spacing. 
    %
    % v = MapGridToBigVector(f2D) -- Maps a quantity defined on the 2D grid
    %                                into vector form.
    %
    % f2D = MapBigVectorToGrid(v) -- Maps a quantity in vector form onto
    %                                the 2D grid. It is often easier to
    %                                acces the desired points in the 2D
    %                                representation.
    %
    % fPara = GetParallelCut(v)   -- Returns the cut in the (positive) 
    %                                parallel direction of a function in
    %                                vector form
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties 
        nP
        nXi
        pMax
        pGridMode
        xiGridMode
        pGridParameter = 0.2 
        pStepParams
        
        stencil = 5;            %3, 5, or 7
        quadrature = 'simpson'; %'simpson' or 'trapezoid'
    end
    
    %%% Grid properties and interface methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    properties 
            
        matSize
        
        %%%%%%%%%%%%%%%%%%
        
        %p grid
        p
        gamma
        dp        
        pWeights
        ddp
        d2dp2
        p2D
        pBig
        gammaBig        
        idsPMax
        
        %%%%%%%%%%%%%%%%%%
        
        %xi grid
        xi
        xiWeights
        ddxi
        d2dxi2
        xi2D
        xiBig
        xi0Id
        
        %%%%%%%%%%%%%%%%%%
             
        %p_para, p_perp space
        pPara2D
        pPerp2D
        
        %%%%%%%%%%%%%%%%%%
        
        %Derivative and integral operators
        ddpMat
        ddxiMat
        d2dp2Mat
        d2dxi2Mat
        d2dpdxiMat
        intdp
        intdxi
        intdpdxi         
        
    end
    
    methods 
        function o = Grid(varargin)
            % Constructor.
            %
            % Usage:
            %   Grid()       -- Initializes an empty object
            %   Gird(oNORSE) -- Uses the settings in the NORSE object
            %
            % Also create a grid with the specified parameters:
            %   Grid(nP,nXi,pMax,pGridMode,xiGridMode)
            %   Grid(nP,nXi,pMax,pGridMode,xiGridMode,pGridParameter)            
            %                        
            
            %Initialize a struct with grid parameters
            o.pStepParams = struct('stepLocation',5,'bulkSpacing',1e-2,...
                                   'Theta',0.01,'stepSharpness',0.1);
            switch nargin
                case 0
                    %An empty object, nothing to do
                case 1
                    if isa(varargin{1},'NORSE')
                        oN = varargin{1};
                        if any([isempty(oN.nP),isempty(oN.nXi),...
                                isempty(oN.pMax),isempty(oN.pGridMode),...
                                isempty(oN.xiGridMode),isempty(oN.pGridParameter)])
                            error('Cannot initialize the grid. Not all the grid parameters have been set.');
                        end
                        o.nP = oN.nP;
                        o.nXi = oN.nXi;
                        o.pMax = oN.pMax;
                        o.pGridMode = oN.pGridMode;
                        o.xiGridMode = oN.xiGridMode;
                        o.pGridParameter = oN.pGridParameter;
                        o.InitializeGrid();
                    else
                       error('The argument must be a NORSE object.'); 
                   end
                case {5,6}
                    o.nP = varargin{1};
                    o.nXi = varargin{2};
                    o.pMax = varargin{3};
                    o.pGridMode = varargin{4};
                    o.xiGridMode = varargin{5};
                    if nargin == 6
                        o.pGridParameter = varargin{6};
                    end
                    o.InitializeGrid();
                otherwise
                    error('Wrong number of input arguments.');
            end            
        end
        
        InitializeGrid(o)
            % Function that constructs the grid using the parameters and
            % settings set as object properties, and calculates all the
            % associated quantities.
        v = MapGridToBigVector(o,f2D)            
            % Reshapes a function given on the 2D grid to a vector
            % representation.
        f2D = MapBigVectorToGrid(o,v)
            % Reshapes a function in vector representation onto the 2D grid.
        fPara = GetParallelCut(o,v)
            % Returns the cut in the (positive) parallel direction of a
            % function in vector form
        Visualize(o)
            % Plots the distribution of the grid points in (p,xi) and
            % (p_para,p_perp) space, as well as the grid-point spacing.
    end

    %%% Internal methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods 
        [x,w,ddx,d2dx2] = CreateUniformGrid(o,n,xMin,xMax)
            % Generates a uniform grid, together with quadrature weights
            % and the first and second finite-diffrences.
        [p,dpds,d2pds2] = DefineStepMapping(o,s)
            % Defines a non-uniform grid mapping with a tanh step in
            % spacing. With this mapping, the spacing is kept small for low
            % p to accurately resolve the bulk, but can be much larger for
            % large p to reduce the computational expense. 
        ConstructBigDifferentiationMatrices(o)
            % Construct differentiation matrices and quadrature weights for
            % the whole grid in vector form.
        m = BuildBigMatrix(o,fOfXi,fOfP,doP0Dependence)
            % Construct a matrix that can be used for operating on a
            % function on the grid (in its vector representation). Gives a
            % vector when multiplying a function on the grid.            
        v = BuildBigVector(o,fOfXi,fOfP)
            % Construct a vector that describes a function on the grid (in
            % its vector representation). Gives a scalar when multiplied by
            % a function on the grid.
        GetSize(o) 
            % Calculates the size in memory of the GRID object by looping
            % through the fields and summing up the variable sizes, and
            % prints it to the console. 
    end
    
    methods (Static) 
        [c,h] = FindOptimalSigmoidMapping(c,d,h,nXi)
            %Finds a grid mapping with a point at xi=0 and endpoints at
            %xi=-1,1, given an initial guess for the parameters c,d and h
            %in the sigmoid mapping for the xi grid. Returns the
            %corresponding values for c and h.                    
        r = Residual(v,d,nXi,id)
            %Returns the residual for the optimization to find a sigmoid
            %grid mapping with a point at xi=0.            
    end      
end
