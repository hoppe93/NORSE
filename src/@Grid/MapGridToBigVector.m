function v = MapGridToBigVector(o,f2D)
    % Reshapes a function given on the 2D grid to a vector
    % representation.
    %
    % Usage:
    %   fV = MapGridToBigVector(f2D)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    v = zeros(o.matSize,1);
    v(1:end-1) = reshape(f2D(2:end,:),o.matSize-1,1); %All points p>0

    %Handle the special point at p=0
    v(end) = f2D(1,o.xi0Id);
        %Use the value at p=0,xi=0. If the function is not single
        %valued at p=0 (this is for instance the case for the xi2D
        %matrix, since it is constructed by concatenating nP copies
        %of the row vector xi), this is the best choice considering
        %the boundary conditions. Otherwise it doesn't matter which
        %xi is used, so we might as well use xi=0.
end
