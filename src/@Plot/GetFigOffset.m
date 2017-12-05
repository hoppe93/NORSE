function figOffset = GetFigOffset(o)
    % Returns the figure offset specified in the NORSE object.
    %
    % Usage:
    %   figOffset = GetFigOffset()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isempty(o.norse)
        figOffset = 0;
    else
        figOffset = o.norse.figOffset; 
    end
end
