function varargout = Dist1D(o,varargin)
    % Plot of the distribution in the parallel direction. Dotted
    % lines denote a negative value of f.
    % 
    % Usage:
    %   PlotDist1D()          -- Stand-alone figure, final time step
    %   PlotDist1D(iteration) -- Stand-alone figure, specified time
    %                            step. Negative values specify
    %                            backwards indexing (-1 gives the
    %                            final time step, -2 the next-to-
    %                            last, and so on).
    %   PlotDist1D(iteration,argsToPlot) -- Assumes this function is
    %           called inside a script which handles all the visual
    %           tweaking. argsToPlot are passed directly to the
    %           plot command and can for instance be used to define
    %           the line style and color. 
    %   h = PlotDist1D(iteration,argsToPlot) -- returns a handle to
    %                                           the created line.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ylims = [1e-20,5];

    %%%%%%%

    oN = o.norse;
    nTimeSteps = numel(oN.times);

    if nargin == 1
        %We're using the routine for a stand-alone plot
        doBarePlot = 0;
        iteration = nTimeSteps; %Use final time step
    elseif nargin == 2                                
        %We're using the routine for a stand-alone plot
        doBarePlot = 0;
        iteration = varargin{1}; %Use a specified time step
    elseif nargin >= 3
        %We're using the routine inside some other script which
        %handles all the visual tweaking
        doBarePlot = 1;   
        iteration = varargin{1};                               
    end
    if iteration <= 0
        iteration = nTimeSteps+1+iteration;
    elseif iteration > nTimeSteps
        warning('The supplied time step is too large. Using the last available time step instead.');
        iteration = nTimeSteps;
    end

    allP = [-flipud(o.grid.p);o.grid.p];            
    fOnGrid = o.grid.MapBigVectorToGrid(oN.f(:,iteration));
    f = [flipud(fOnGrid(:,1));fOnGrid(:,end)]; %f at xi=-1 and xi=1
    if doBarePlot 
        h = semilogy(allP,f,varargin{2:end});                
    else
        %Prepare the inital state and the initial Maxwellian
        fInitial= o.grid.MapBigVectorToGrid(oN.f(:,1));
        fInitial = [flipud(fInitial(:,1));fInitial(:,end)];
        fMaxw = oN.maxwellianPreFactor(0) * ...
                              exp( (1-o.grid.gamma)/(oN.Theta(0)) );
        fMaxw = [flipud(fMaxw);fMaxw];

        %Initialize the figure
        figure(o.GetFigId('Dist1D'));
        clf;

        %Plot
        semilogy(allP,fInitial,'--k');
        hold on;
        semilogy(allP,fMaxw,'-.r');      
        semilogy(oN.runawayRegion.pcs(end)*[1,1],ylims,':m');                
        semilogy(allP,f,'-b');
        semilogy(allP,-f,':b');

        %Label the plot
        hL = legend('Initial','Maxwellian at $t_0$','$p_c(\xi=1)$');
        hL.Interpreter = 'latex';
        xlabel('$p_{||}$','Interpreter','latex')
        ylabel('$F$ ($p_{\perp}=0$)','Interpreter','latex')
        ylim(ylims);
        title(sprintf('$\\tau = %.3g$ (relativistic coll. times)',...
                        oN.times(iteration)),'Interpreter','latex');
    end

    if nargout == 1
        varargout{1} = h;
    else
        varargout = {};
    end
end
