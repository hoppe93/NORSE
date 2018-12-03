function PerformCalculation(o, varargin)
    %Carries out the various steps in a standard NORSE calculation. Optionally 
    %uses an existing Grid object set in the NORSE class property grid -- 
    %otherwise a new grid is generated using the settings in the NORSE object.
    %
    % Usage:
    %   PerformCalculation() 
    %   PerformCalculation(useExistingGrid)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    tStart = tic;
    %separate cases according to the need for external distribuiton
    switch o.initialDistribution
        case 0
            if nargin==2 && varargin{1}
                %Use an existing grid object, but make sure it is initialized
                o.grid.InitializeGrid();
            else
                %Create a new grid object
                o.grid = Grid(o);
            end
            o.Initialize();
            o.timeAdvance.AdvanceInTime();                    
            o.PostProcess();
            o.timing.total = toc(tStart);
            o.PrintTimings();
        case 1
            if nargin==2 && varargin{1}
                %Use an existing grid object, but make sure it is initialized
                o.grid.InitializeGrid();
            else
                %Create a new grid object
                o.grid = Grid(o);
            end
            o.Initialize();
            o.timeAdvance.AdvanceInTime();                    
            o.PostProcess();
            o.timing.total = toc(tStart);
            o.PrintTimings();
        case 2
            if nargin==2 && varargin{1}
                %Use an existing grid object, but make sure it is initialized
                o.grid.InitializeGrid();
            else
                %Create a new grid object
                o.grid = Grid(o);
            end
            o.Initialize();
            o.timeAdvance.AdvanceInTime();                    
            o.PostProcess();
            o.timing.total = toc(tStart);
            o.PrintTimings();
        case 3
            if nargin==2 && varargin{1}
                %Use an existing grid object, but make sure it is initialized
                o.grid.InitializeGrid();
            else
                %Create a new grid object
                o.grid = Grid(o);
            end
            o.Initialize();
            o.timeAdvance.AdvanceInTime();                    
            o.PostProcess();
            o.timing.total = toc(tStart);
            o.PrintTimings();
        case 4
            if nargin==3 && varargin{1}
                %Use an existing grid object, but make sure it is initialized
                o.grid.InitializeGrid();
                o.Initialize(varargin{2});
            elseif nargin==2 
                %Create a new grid object
                o.grid = Grid(o);
                o.Initialize(varargin{1});
            end
            o.timeAdvance.AdvanceInTime();                    
            o.PostProcess();
            o.timing.total = toc(tStart);
            o.PrintTimings();
        otherwise
            error('Invalid initial distribution.');
    end
end
