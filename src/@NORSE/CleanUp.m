function CleanUp(o)
    % Clear large internal variables that are not needed when
    % saving the NORSE object. The object can still be used to
    % continue the calculation (after calling 
    % RebuildInternalVariables).
    %
    % Usage:
    %   CleanUp()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    o.potentials        = [];
    o.collisionOperator = [];
    o.heatSink          = [];
    o.particleSource    = [];

    o.EFieldOperator = [];            
    o.synchrotronOperator = [];
    
    o.timeAdvance.CleanUp();
    
    o.Upsilon0 = [];
    o.Upsilon1 = [];
    o.Upsilon2 = [];
    o.Pi       = [];
end 
