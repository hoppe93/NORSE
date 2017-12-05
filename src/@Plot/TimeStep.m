function TimeStep(o)
    % Plots the time step used by the adaptive time step scheme, as
    % a function of time (normalized to the initial time step).
    % Also plots the number of GMRES iterations used to determine
    % the time step evolution.
    %
    % Usage:
    %   TimeStep()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if o.norse.timeAdvanceMode == 2
        figure(o.GetFigId('TimeStep'));
        clf;                                

        % Plot time steps used
        subplot(2,1,1)
        t = o.norse.timeAdvance.allTimes;
        plot(t,o.norse.timeAdvance.dtsUsed/o.norse.dt,'.-','ButtonDownFcn',@o.PlotAtTimeX);
        ylabel('$\mathrm{d}\tau/\mathrm{d}\tau_0$','Interpreter','latex');
        title(sprintf('$\\mathrm{d}\\tau_0 = %.2g \\tau$',o.norse.dt),...
                                            'Interpreter','latex');

        % Plot number of GMRES iterations
        subplot(2,1,2)
        plot(t,o.norse.timeAdvance.gmresIters,'.-')
        ylabel('GMRES iterations','Interpreter','latex')
        xlabel('$\tau$','Interpreter','latex');
    else
        fprintf('Nothing to show, since a constant time step was used.\n'); 
    end
end
