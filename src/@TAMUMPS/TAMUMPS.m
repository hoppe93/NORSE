classdef TAMUMPS < TimeAdvance
    % TAMUMPS -- Uses the MUMPS library (direct solver) to advance the
    %            NORSE system in time.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Usage:
    % %%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Initialization:
    %       TAMUMPS(mumpsPath) -- empty object
    %
    %   In NORSE:
    %       TAMUMPS(mumpsPath,oNORSE)
    %       AdvanceInTime()
    %
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Interface methods and output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods
        AdvanceInTime(o,varargin)
            % Advance the system in time using the MUMPS solver.
    end

    %%% Internal properties and methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        mumpsInfo

        % %%%%%%%%%%%%%%%
        % MUMPS options
        % %%%%%%%%%%%%%%%

        precinv
        % Stopping criterion for iterative refinement
        %   Possible values: >= 0
        %                    values < 0 translates to sqrt(eps)
        %   Default value:   sqrt(eps)
        % (more info: CNTL(2) in MUMPS manual)

        thresholdPivoting
        % Relative threshold for numerical pivoting
        %   Possible values:  < 0  let MUMPS decide!
        %                     = 0  no pivoting performed
        %                     > 0  numerical pivoting performed
        %                          For unsymmetric matrices, values > 1 are treated as 1
        %                          For symmetric matrices, values > 0.5 are treated as 0.5
        %   Default value:   1e-2  For unsymmetric or general symmetric matrices
        %                       0  For symmetric positive definite matrices
    end

    methods
        function o = TAMUMPS(mumpsPath,varargin)
            % Constructor.
            %
            % Usage:
            %   TAMUMPS()
            %   TAMUMPS(mumpsPath,oNORSE) -- call to use in NORSE.
            %                                Specifies also path to MUMPS
            %                                MEX file (may be empty).

            o@TimeAdvance(varargin{:});
            o.CalculateNumberOfSteps();

            % Initialize MUMPS
            o.InitializeMUMPS(mumpsPath);

            % Set default properties
            % < 0 means sqrt(eps)
            o.precinv = -1;

            % < 0 means MUMPS will get to decide automatically
            o.thresholdPivoting = -1;

        end
        
        InitializeMUMPS(o,mumpsPath);
        % Initialize the MUMPS solver
        
        MUMPSAnalyze(o,matrix);
        % Run the MUMPS analysis stage on the matrix (required before
        % solve)
        
        f = MUMPSSolve(o,matrix,rhs);
        % Solve the linear system using MUMPS
        
    end

end
