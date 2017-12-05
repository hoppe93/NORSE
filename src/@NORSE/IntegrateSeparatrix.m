function [xiS,pS] = IntegrateSeparatrix(o,pc,iteration)
    % Strating from p=pc, integrate 'backwards' in xi (from xi=1 to
    % xi=-1) to trace out the entire separatrix. The returned xis
    % are points on the xi grid, but are in reverse
    % order, and not all points are necessarily included. iteration
    % specifies the time step to use for the collisional friction.
    % The vectors xis and ps together define a curve that is the
    % separatrix.
    %
    % Usage:
    %   [xis,ps] = IntegrateSeparatrix(pc,iteration)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    time = o.times(iteration);            
    t0 = o.referenceTime;

    %Make a spline of the potential Pi and take its derivatives.
    Pi = o.grid.MapBigVectorToGrid(o.Pi(:,iteration));
    PiSpline2D = csapi({o.grid.p,o.grid.xi},Pi);
    dPidp = fndir(PiSpline2D,[1;0]);
    dPidxi = fndir(PiSpline2D,[0;1]);            

    %Set up the integrand and integration parameters
    TNorm = o.lnLambdaBar(t0)/(o.Theta(t0)*o.kappa(t0));
    integrand = @(xi,p) o.SeparatrixIntegrationEquation(xi,p,TNorm,...
                            o.EHat(time),o.sigma(time),dPidp,dPidxi);
    interruptIntegrationFcn = @(xi,p,flag) o.InterruptIntegration(...
                                                p,flag,o.grid.pMax);
    opts = odeset('OutputFcn',interruptIntegrationFcn,'RelTol',1e-5);                

    %%% Integrate the trajectory %%%
    %Turn off a common warning associated with the separatrix
    %"leaving" the numerical grid
    warning off MATLAB:ode45:IntegrationTolNotMet
    xiSpan = flipud(o.grid.xi(2:end));
    xiSpan(1) = 1-1e-6; 
    [xiS,pS] = ode45(integrand,xiSpan,pc,opts); 
                            %This is slow. The main cause is the
                            %evaluation of the 2D spline.
    warning on MATLAB:ode45:IntegrationTolNotMet

    %We didn't quite start the integration at xi=1, so overwrite
    %the first value
    xiS(1) = 1;
    pS(1) = pc;
end
