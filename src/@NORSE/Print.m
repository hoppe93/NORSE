function Print(o,str,varargin)
    %Prints text to the console if silent mode is not
    %enabled. Uses the syntax of fprintf.
    %
    % Usage:
    %   Print(str)
    %   Print(str,argsToFprintf)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~o.silent
        fprintf(str,varargin{:});
    end
end
