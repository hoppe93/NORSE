function fPara = GetParallelCut(o,v)
    % Returns the cut in the (positive) parallel direction of a function
    % in vector form
    %
    % Usage:
    %   fPara = GetParallelCut(fV)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    f2D   = o.MapBigVectorToGrid(v);
    fPara = f2D(:,end);
end
