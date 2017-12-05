function DistMovie1D(o,varargin)
    % Plays a "movie" of the distribution evolution in the parallel
    % direction. yLimits is a vector specifying the axis limits for
    % the distribution F (on the format used by ylim()) and tTot is
    % the duration of the movie in second (determines the refresh
    % rate). Default: yLimits = [1e-20,5], tTot = 5 (s)
    %
    % Note that the movie plots the saved time steps with equal
    % refresh rate. Depending on how the save steps are chose, the
    % movie might not be linear in physical time (seconds).
    %
    % Usage:
    %   DistMovie1D()
    %   DistMovie1D(yLimits)
    %   DistMovie1D(yLimits,tTot)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    yL   = [1e-20,5]; %Plot limits
    tTot = 5;         %Total time (s) of movie            

    %%%%%%%%%%

    %Handle inputs
    switch nargin
        case 1
            %nothing to do
        case 2
            yL = varargin{1};
        case 3
            yL = varargin{1};
            tTot = varargin{2};
        otherwise
            error('Wrong number of input arguments');
    end

    oN = o.norse;
    nSteps = numel(oN.times);
    pauseTime = tTot/nSteps;

    %Prepare grid, initial and Maxwellian distributions
    allP = [-flipud(o.grid.p);o.grid.p];                        
    fInitial= o.grid.MapBigVectorToGrid(oN.f(:,1));
    fI = [flipud(fInitial(:,1));fInitial(:,end)]; %combine xi=-1 and xi=1
    fMaxw = exp( (1-o.grid.gamma)/(oN.Theta(oN.referenceTime)) );
    fM = [flipud(fMaxw);fMaxw];

    %Prepare figure
    figure(o.GetFigId('DistMovie1D'));            
    clf;
    hA = axes();

    %Plot initial state
    semilogy(allP,fI,'--k');
    hold on;
    semilogy(allP,fM,'-.r');          
    hPos = semilogy(allP,fI,'-b');
    hNeg = semilogy(allP,-fI,':b');

    %Label the plot            
    hL = legend('Initial','Maxwellian at $t_0$');
    hL.Interpreter = 'latex';
    xlabel('$p_{||}$','Interpreter','latex')
    ylabel('$F$ ($p_{\perp}=0$)','Interpreter','latex')
    ylim(yL);            
    hTitle = title(sprintf('$\\tau = %.2g$ (relativistic coll. times)',...
                                oN.times(1)),'Interpreter','latex');

    %Set the lower x limit based on the extend of the dist in the
    %xi=-1 direction
    idXMin = min([find(fI>yL(1),1),o.grid.nP]);
    xlim([allP(idXMin),allP(end)]);

    %Plot each time step in turn
    for i = 2:nSteps
        %Update the dist
        fOnGrid = o.grid.MapBigVectorToGrid(oN.f(:,i));
        f = [flipud(fOnGrid(:,1));fOnGrid(:,end)];
        hPos.YData = f;
        hNeg.YData = -f;                

        %Update the lower x limit
        idXMin = min([find(f>yL(1),1),o.grid.nP]);
        if allP(idXMin) < hA.XLim(1)
            hA.XLim(1) = allP(idXMin);
        end

        hTitle.String = sprintf('$\\tau = %.3g$ (relativistic coll. times)'...
                                                        ,oN.times(i));
        drawnow;
        pause(pauseTime);
    end            
end
