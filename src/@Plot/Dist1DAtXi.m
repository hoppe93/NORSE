function Dist1DAtXi(o,idT,idXi,varargin)
    % Plots the distribution at a specified xi (id on the xi-grid)
    % at a given time step. If idTime=-1, the final time step is
    % used, if idT=-2, the next-to-final, and so on. The optional
    % argument yLimits is a vector specifying the plot limits on
    % the value of f. Default: yLimits = [1e-20,5]
    %
    % Usage:
    %   Dist1DAtXi(idTime,idXi)
    %   Dist1DAtXi(idTime,idXi,yLimits)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Defaults
    yL = [1e-20,5]; 

    %%%%%%%%%%%

    %Handle input
    oN = o.norse;
    switch nargin
        case 3
            %Nothing to do
        case 4
            yL = varargin{1};
        otherwise
            error('Invalid number of input arguments.');
    end
    if idT < 0
        idT = numel(oN.times) + 1 + idT;
    end

    %Prepare distribution

    p = o.grid.p;
    fOnGrid = o.grid.MapBigVectorToGrid(oN.f(:,idT));
    fInitial = o.grid.MapBigVectorToGrid(oN.f(:,1));
    f = fOnGrid(:,idXi);
    fI = fInitial(:,idXi);

    %Prepare figure
    figure(o.GetFigId('Dist1DAtXi'));
    clf;

    %Plot
    semilogy(p,fI,'--k');            
    hold on;
    semilogy(p,f,'-b');
    semilogy(p,-f,':b');

    %Label the plot
    hF = legend('Initial','$F$');
    hF.Interpreter = 'latex';
    xlabel('$p$','Interpreter','latex');
    ylabel('$F$','Interpreter','latex');
    ylim(yL);
    title(sprintf('$\\tau = %.2g$ (relativistic coll. times), $\\xi = %.4g$',...
                    oN.times(idT),o.grid.xi(idXi)),'Interpreter','latex');
end
