function oCopy = copyElement(o)
    % Overrides Matlab's standard functionality for making
    % (shallow) copies of handle objects to ensure that copies of
    % the referenced handle objects are also made.
    %
    % Usage:
    %      oCopy = copy(o)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Make a shallow copy using the standard copy functionality
    oCopy = copyElement@matlab.mixin.Copyable(o);

    %Determine what handle objects are instantiated
    isGrid   = ~isempty(o.grid);
    isPot    = ~isempty(o.potentials);
    isCollOp = ~isempty(o.collisionOperator);
    isHS     = ~isempty(o.heatSink);
    isPS     = ~isempty(o.particleSource);
    isTA     = ~isempty(o.timeAdvance);
    isPlot   = ~isempty(o.plot);

    %Make deep copies of the handle objects
    if isGrid
        oCopy.grid              = copy(o.grid);
    end
    if isPot
        oCopy.potentials        = copy(o.potentials);
    end
    if isCollOp
        oCopy.collisionOperator = copy(o.collisionOperator);
    end
    if isHS
        oCopy.heatSink          = copy(o.heatSink);
    end
    if isPS
        oCopy.particleSource    = copy(o.particleSource);                          
    end
    if isTA
        oCopy.timeAdvance       = copy(o.timeAdvance); 
    end
    if isPlot
        oCopy.plot              = copy(o.plot);             
    end

    %Also update the references between the different objects
    if isGrid
        if isPot
            oCopy.potentials.norse        = oCopy;
            oCopy.potentials.grid         = oCopy.grid;                
        end
        if isCollOp && isPot
            oCopy.collisionOperator.norse = oCopy;
            oCopy.collisionOperator.grid  = oCopy.grid;
            oCopy.collisionOperator.potentials = oCopy.potentials;
        end
        if isHS && isPS
            oCopy.heatSink.norse          = oCopy;
            oCopy.heatSink.grid           = oCopy.grid;
            oCopy.heatSink.particleSource = oCopy.particleSource;
        end
        if isPS
            oCopy.particleSource.norse    = oCopy;
            oCopy.particleSource.grid     = oCopy.grid;
        end 
        if isTA
            oCopy.timeAdvance.norse       = oCopy;            
        end
        if isPlot
            oCopy.plot.norse              = oCopy;
            oCopy.plot.grid               = oCopy.grid;
        end
    end
end 
