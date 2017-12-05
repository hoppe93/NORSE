function [p,dpds,d2pds2] = DefineStepMapping(o,s)
    % Defines a non-uniform grid mapping with a tanh step in
    % spacing. With this mapping, the spacing is kept small for low
    % p to accurately resolve the bulk, but can be much larger for
    % large p to reduce the computational expense. 
    %
    % The grid spacing has the form
    %   dpds = a*tanh((s-sp)/kappa) + b,
    % which can be integrated to get the expression for p(s). The
    % exact shape of the resulting grid is governed by the settings
    % in the struct pStepParams. kappa is given by the struct field
    % stepSharpness, whereas sp is calculated from a combination of
    % stepSharpness, stepLocation, bulkSpacing and Theta. a and b
    % are determined by solving a system of constraint equations,
    % given kappa, sp, s, nP, pMax and stepLocation.
    %
    % Usage:
    %   [p,dpds,d2pds2] = DefineStepMapping(s)
    %
    % s should be a uniform grid on [0,1].
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Obtain settings from the pStepParams struct                        
    yToP = sqrt(2*o.pStepParams.Theta); %Conversion between y and p
    pStep = yToP*o.pStepParams.stepLocation; 
    deltaP0 = yToP*o.pStepParams.bulkSpacing;
    kappa = o.pStepParams.stepSharpness; 

    % Determine derived quantities
    maxFractionOfPoints = 0.6; %to use for the bulk
    N0 = round(min(pStep/deltaP0,maxFractionOfPoints*o.nP)); %points in bulk
    Np = round(N0 + 2.5*kappa*o.nP); %points in bulk and step regions
    sp = s(Np);                        
    sN0 = s(N0);

    if sN0 == 0
        error('Unsuitable grid parameters.');
    end

    %%% Solve for suitable grid parameters a and b 
    % These expressions follow from the constraints 
    %   p(s=1) = pMax
    %   p(N0)  = pStep
    K11 = kappa * ( log(cosh((1-sp)/kappa)) - log(cosh(sp/kappa)) );
    K12 = 1;
    K21 = kappa * ( log(cosh((sN0-sp)/kappa)) - log(cosh(sp/kappa)) );
    K22 = sN0;
    K = [K11,K12;K21,K22];
    if any(isinf(K(:))) || any(isnan(K(:)))
        error('Unsuitable grid parameters.');
    end

    ab = K\[o.pMax;pStep];
    a = ab(1);
    b = ab(2);                        
    if b<a
        %This is required to guarantee positivity of dpds (tanh ->
        %-1 as the argument -> 0)
        b=a+eps;
    end

    %%% Define the grid                                                
    p = a*kappa*( log(cosh((s-sp)/kappa)) - log(cosh((sp)/kappa)) ) + b*s;
    dpds = a*tanh((s-sp)/kappa) + b;
    d2pds2 = (a/kappa)./cosh((s-sp)/kappa).^2;

    if (p(end)-p(end-1)) < (p(2)-p(1))
        error('The grid will not have the desired properties for the chosen settings. Decrease nP or bulkSpacing!' );
    end
    if any(dpds<=0)
        error('The grid-point slope is not strictly positive! The p grid mapping parameters are invalid.');
    end                   
end
