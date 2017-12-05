function v = BuildBigVector(o,fOfXi,fOfP)
    % Construct a vector that describes a function on the grid (in
    % its vector representation). Gives a scalar when multiplied by
    % a function on the grid.
    %
    % Usage:
    %   v = BuildBigVector(fOfXi,fOfP)
    %
    % fOfXi describes the xi dependence and has dimension [1,nXi].
    % Similarly, fOfP describes the p dependence and is [1,nP].
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    v = zeros(o.matSize,1);
    v(1:(end-1)) = kron(fOfXi,fOfP(2:end)); %all points with p>0
    v(end) = fOfP(1); %At p=0. This cannot depend on xi. Is it then OK to only use the p value? At least it is in the intVector case since fOfP=0 at p=0
end 
