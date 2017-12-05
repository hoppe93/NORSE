function o = Rescale(o,tFac)
    % Rescales the time axis so that tNew = tFac*tOld.
    %
    % Usage:
    %   par = par.Rescale(tFac)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~o.isScalar
        %Change the times at which the spline is defined
        breaks = tFac*o.parameterSpline.breaks;
        o.parameterSpline.breaks = breaks;
        o.tMin = breaks(1);
        o.tMax = breaks(end);

        %Change the polynomial coefficients to match the new times
        polynomialOrder = o.parameterSpline.order-1;
        scaleFacs = tFac.^(polynomialOrder:-1:0);
        o.parameterSpline.coefs = o.parameterSpline.coefs*diag(1./scaleFacs);                
    end            
end
