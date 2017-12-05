function GMRESInfo(o)
    % Plots information about the performance of GMRES in the
    % iterative and adaptive time-advancement schemes. Data is only
    % available for the time steps when GMRES was used.
    %
    % Flag:       Error flag (0: OK)
    % Residual:   Value of residual at convergence (or abortion)
    % Iterations: Number of iterations to convergence
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    oN  = o.norse;
    oTA = oN.timeAdvance;
    if ~oN.timeAdvanceMode
       warning('The time-evolution was performed using the direct solver -- no information to plot.');
       return
    end

    figure(o.GetFigId('GMRESInfo'));
    clf;

    subplot(3,1,1)
    plot(oTA.gmresFlags,'.-')
    ylabel('GMRES flag','Interpreter','latex')
    title('1: Did not converge, 2: Preconditioner illconditioned, 3: Stagnated',...
                                            'Interpreter','latex');

    subplot(3,1,2)
    semilogy(oTA.gmresRess,'.-')
    ylabel('GMRES residual','Interpreter','latex')

    subplot(3,1,3)
    plot(oTA.gmresIters,'.-')
    ylabel('GMRES iterations','Interpreter','latex')
    xlabel('step','Interpreter','latex')
    drawnow
end
