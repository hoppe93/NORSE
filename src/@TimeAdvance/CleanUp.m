function CleanUp(o)
    % Clear large internal variables that are not needed when
    % saving the NORSE object.
    % Usage:
    %   CleanUp()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    if isvalid(o) %If the object has not been removed, clear some fields
        o.matrix = [];    
        o.rhs = [];
    end
        
