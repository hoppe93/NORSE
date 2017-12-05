function DistMovie2D(o,varargin)
    % Plays a "movie" of the 2D distribution evolution in
    % (p_para,p_perp) space. 
    %
    % Note that the movie plots the saved time steps with equal
    % refresh rate. Depending on how the save steps are chose, the
    % movie might not be linear in physical time (seconds).
    %
    % Usage:
    %   DistMovie2D()
    %   DistMovie2D(contours)            
    %   DistMovie2D(contours,autoTrim)
    %
    % contours is a vector specifying the contour values to use in
    % contourf (a logarithmic scale), and autoTrim specifies
    % whether to automatically crop the plots to remove excess
    % white space in the p_perp and negativ p_para directions.
    % Default: contours = linspace(-20,1,15), autoTrim = true
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Default settings
    contours    = linspace(-20,1,15); %Contours values to plot            
    useAutoTrim = true;
    pauseTime   = 0; %The movie can be slowed using a positive value
    mirrorDist  = 0; %Mirrors the plot in the parallel axis
    skip        = 1; %Skips time steps to speed up the movie

    %%%%%%%%%%

    %Handle inputs
    switch nargin
        case 1
            %nothing to do
        case 2
            contours = varargin{1};
        case 3
            contours = varargin{1};
            useAutoTrim = varargin{2};                
        otherwise
            error('Wrong number of input arguments');
    end

    oN        = o.norse;
    nSteps    = numel(oN.times);            
    pMax      = o.grid.pMax;
    cutOff    = contours(1);

    hF        = figure(o.GetFigId('DistMovie2D'));
    colormap(o.ColorMap(numel(contours)));
    pCutNeg   = 0;
    pPerpCut  = 0;  

    %Make the grid
    paraGrid = o.grid.pPara2D;
    perpGrid = o.grid.pPerp2D;
    if mirrorDist
        paraGrid = [flipud(paraGrid);paraGrid];
        perpGrid = [-flipud(perpGrid);perpGrid];
    end


    for i = [1:skip:(nSteps-1),nSteps]
        clf(hF);        

        %Get the dist
        f      = o.grid.MapBigVectorToGrid(oN.f(:,i));
        f(f<0) = realmin; %disregard negative values for log plotting
        f      = log10(f);
        if mirrorDist
            f = [flipud(f);f];
        end

        %Plot contours and label the plot 
        contourf(paraGrid,perpGrid,f,contours);
        xlabel('$p_{||}$','Interpreter','latex')
        ylabel('$p_{\perp}$','Interpreter','latex')
        title(sprintf('$\\tau=%.3f$ (relativistic coll. times)',...
                                oN.times(i)),'Interpreter','latex');
        hC = colorbar;
        hC.Label.String = '$\log_{10}(F)$';
        hC.Label.Interpreter = 'Latex';
        axis equal

        %Make sure the colour scale is consistent
        hC.Limits = [contours(1),contours(end)];
        caxis([contours(1),contours(end)]);

        %Determine the plot limits
        if useAutoTrim 
            %--Negative xi
            pCutNeg  = min([pCutNeg,1.05*min(paraGrid(f>cutOff))]);
            %--p_perp
            pPerpCut = max([pPerpCut,1.05*max(perpGrid(f>cutOff))]);
        else
            pCutNeg  = -pMax;
            pPerpCut = pMax;
        end
        xlim([pCutNeg,pMax]);
        ylim([0,pPerpCut]);

        drawnow
        pause(pauseTime);
    end
end        
