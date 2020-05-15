function ConstructBigDifferentiationMatrices(o)
    % Construct differentiation matrices and quadrature weights for
    % the whole grid in vector form.
    %
    % Usage:
    %   ConstructBigDifferentiationMatrices()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Differentiation matrices
    o.ddpMat = o.BuildBigMatrix(speye(o.nXi),o.ddp,1);
    o.ddxiMat = o.BuildBigMatrix(o.ddxi,speye(o.nP),0);
    o.d2dp2Mat = o.BuildBigMatrix(speye(o.nXi),o.d2dp2,1);
    o.d2dxi2Mat = o.BuildBigMatrix(o.d2dxi2,speye(o.nP),0);
    o.d2dpdxiMat = o.BuildBigMatrix(o.ddxi,o.ddp,0);    

    %Integrals over the grid (including Jacobian)   
    o.intdp = o.BuildBigVector(ones(size(o.xi))',4*pi*o.p'.^2.*o.pWeights');
    o.intdpdxi = o.BuildBigVector(o.xiWeights',2*pi*o.p'.^2.*o.pWeights');
end 
