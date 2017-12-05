function f2D = MapBigVectorToGrid(o,v)
    % Reshapes a function in vector representation onto the 2D grid.
    %
    % Usage:
    %   f2D = MapBigVectorToGrid(fV)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    f2D = zeros(o.nP,o.nXi);
    f2D(2:end,:) = reshape(v(1:end-1),o.nP-1,o.nXi); %All points p>0
    f2D(1,:) = v(end); %At p=0, the value should be the same for all xi
end
