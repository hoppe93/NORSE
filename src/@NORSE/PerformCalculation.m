function PerformCalculation(o, varargin)
    %Carries out the various steps in a standard NORSE calculation. Optionally 
    %uses an existing Grid object set in the NORSE class property grid -- 
    %otherwise a new grid is generated using the settings in the NORSE object.
    %
    % Usage:
    %   PerformCalculation() 
    %   PerformCalculation(useExistingGrid)
    %   PerformCalculation(useExternalInput)
    %   PerformCalculation(useExistingGrid, useExternalInput)
    %
    % The useExternalInput variable is a structure with fields
    % corresponding to the external distribution and external grid vectors.
    % For more details, please see Initialize.m.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    tStart = tic;
    %Separate cases according to the need for external distribuiton
    if o.initialDistribution >= 0 && o.initialDistribution <=3
        if nargin==2 && varargin{1}
            %Use an existing grid object, but make sure it is initialized
            o.grid.InitializeGrid();
        else
            %Create a new grid object
            o.grid = Grid(o);
        end
        o.Initialize();
    elseif o.initialDistribution == 4
        % This case requires external distribution and grid vectors given
        % to the Initialize method. See Initialize.m for more details.
        if nargin==3 && varargin{1}
            %Use an existing grid object, but make sure it is initialized
            o.grid.InitializeGrid();
        elseif nargin==2 
            %Create a new grid object
            o.grid = Grid(o);
        else
            error('Invalid number of input arguments.');
        end
        % The last input argument is the external data
        o.Initialize(varargin{end});
            
    else
       error('Invalid initial distribution.');
    end
    o.timeAdvance.AdvanceInTime();                    
    o.PostProcess();
    o.timing.total = toc(tStart);
    o.PrintTimings();
end
