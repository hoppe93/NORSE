function Visualize(o)
    % Plots the distribution of the grid points in (p,xi) and
    % (p_para,p_perp) space, as well as the grid-point spacing.
    %
    % Usage:
    %   Visualize()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Figure of spatial distribution of grid points
    hF = figure();
    clf;
    hF.Name = 'Spatial distribution of grid points';

    subplot(2,1,1) % (p,xi)-space
    scatter(o.xi2D(:),o.p2D(:),5,'k','filled');
    xlabel('$\xi$','Interpreter','latex')
    ylabel('$p$','Interpreter','latex')
    title('\bf($p$,$\xi$)-space','Interpreter','latex');

    subplot(2,1,2) %(p_para,p_perp)-space
    scatter(o.pPara2D(:),o.pPerp2D(:),5,'k','filled');
    xlabel('$p_{||}$','Interpreter','latex')
    ylabel('$p_{\perp}$','Interpreter','latex')                        
    title('\bf($p_{||}$, $p_{\perp}$)-space','Interpreter','latex');
    axis equal


    %Figure of grid-point spacing
    hF = figure();
    clf;
    hF.Name = 'Grid-point spacing';

    subplot(2,1,1);
    plot(o.p,'k')
    xlabel('Grid point','Interpreter','latex')
    ylabel('$p$','Interpreter','latex')

    subplot(2,1,2);
    plot(o.xi,'k')
    xlabel('Grid point','Interpreter','latex')
    ylabel('$\xi$','Interpreter','latex')
end
