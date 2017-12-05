function PerformCalculation(o, varargin)
    %Carries out the various steps in a standard NORSE calculation. Optionally 
    %uses and existing Grid object set in the NORSE class property grid -- 
    %otherwise a new grid is generated using the settings in the NORSE object.
    %
    % Usage:
    %   PerformCalculation() 
    %   PerformCalculation(useExistingGrid)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    tStart = tic;
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
end
