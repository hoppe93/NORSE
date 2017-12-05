function DistVsTime(o,varargin)
    % Contour plot of the distribution in the parallel direction as
    % a function of time. 
    %
    % Clicking at a point in the plot produces a new 1D plot with
    % the paralllel distribution at that time.
    % 
    % Usage:
    %   DistVsTime()         
    %   DistVsTime(contours) 
    %   DistVsTime(contours,autoTrim)
    %
    % contours is a vector of contour values of log10(f) to plot
    % and autoTrim specifies whether to automatically crop the plot
    % to not show grid points where f<contours(1) at all times.
    % Default: linspace(-12,fMax,35), autoTrim = true
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Defaults
    cutoff = 1e-12;
    contours = [];
    useAutoTrim = 1;

    %Handle inputs
    switch nargin
        case 1
            %Nothing to do
        case 2
            contours = varargin{1};
            cutoff = 10^contours(1);
        case 3
            contours = varargin{1};
            useAutoTrim = varargin{2};
            cutoff = 10^contours(1);
        otherwise
            error('Invalid number of input arguments.');
    end


    %%% Prepare the parallel distribution at all times %%%
    oN = o.norse;
    allP = [-flipud(o.grid.p);o.grid.p(2:end)];
    fOfT = zeros(2*o.grid.nP-1,oN.nSaveSteps);

    % The f array has the f at different times in different
    % columns. Each row thus gives f at a specific grid point at
    % all times. We just need to find the appropriate rows for the
    % grid points we want. 
    fOfT(1:o.grid.nP-1,:) = flipud(oN.f(1:o.grid.nP-1,:)); %xi=-1
    fOfT(o.grid.nP,:) = oN.f(end,:); %p=0 (all xi)
    fOfT(o.grid.nP+1:end,:) = oN.f((end-(o.grid.nP-1)):(end-1),:); %xi=1
    fOfT(fOfT<=cutoff) = cutoff;

    %Remove grid points for high or low p where f<cutoff at all times
    if useAutoTrim                
        ids = all(fOfT'==cutoff); %Grid points (columns) where f<cutoff at all times
        idPMin = find(ids==0,1)-1; 
        idPMax = find(ids==0,1,'last');     
        if isempty(idPMin) || ~idPMin
            idPMin = 1;
        end
        if isempty(idPMax) || idPMax == size(fOfT,1)
            idPMax = size(fOfT,1);
        else
            idPMax = idPMax+1;
        end

        %Remove the unwanted grid points
        allP = allP(idPMin:idPMax);
        fOfT = fOfT(idPMin:idPMax,:);
    end

    %Calculate contours (if not specified)
    if isempty(contours)
        contours = linspace(log10(cutoff),max(fOfT(:)),35);
    end


    %%% Plot it %%%
    figure(o.GetFigId('DistVsTime'));         
    %Flip the f array, so as to get time on the y axis
    contourf(allP,oN.times,log10(fOfT'),contours,...
             'LineColor','none','ButtonDownFcn',@o.PlotAtTimeY);                        
    xlabel('$p_{||}$','Interpreter','latex');
    ylabel('$\tau$','Interpreter','latex');
    colormap(o.ColorMap(numel(contours)));
    hC = colorbar;
    hC.Label.String = '$\log_{10}(F)$';
    hC.Label.Interpreter = 'Latex';
end
