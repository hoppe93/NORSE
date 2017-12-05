function varargout = GetRelativisticNormalization(o,tMax,dt,T,varargin)
    %Converts tMax & dt given in thermal collision times to the
    %relativistic collision times used in NORSE. Optionally also
    %returns the scale factor used. Optionally converts the grid
    %max in thermal momenta (y=gamma*v/v_th) to the value in p. The
    %temperature should be in eV.
    %
    % Usage: 
    %   [tMax,dt] = GetRelativisticNormalization(tMax,dt,T)
    %   [tMax,dt,tScaleFactor] = GetRelativisticNormalization(tMax,dt,T)
    %   [tMax,dt,pMax] = GetRelativisticNormalization(tMax,dt,T,yMax)
    %   [tMax,dt,pMax,tScaleFactor] = GetRelativisticNormalization(tMax,dt,T,yMax)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Theta = T*(o.constants.e/(o.constants.m*(o.constants.c)^2)); %T in eV

    normFac = (2*Theta).^(3/2);
    tMax = tMax.*normFac;
    dt = dt.*normFac;

    if nargin == 5
        yMax = varargin{1};
        pMax = sqrt(2*Theta).*yMax;
        if nargout == 3
            varargout = {tMax,dt,pMax};
        elseif nargout == 4
            varargout = {tMax,dt,pMax,normFac};
        else
            error('Invalid number of output arguments.');
        end
    elseif nargout == 2                
        varargout = {tMax,dt};
    elseif nargout == 3
        varargout = {tMax,dt,normFac};
    else
        error('Invalid combination of input and output arguments.');
    end

end
