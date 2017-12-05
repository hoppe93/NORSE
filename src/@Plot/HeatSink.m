function HeatSink(o,varargin)
    % Shows the magnitude of the heat sink, as well as the
    % contributions from its various components, as functions of
    % time. Also shows the momentum-space shape of the heat sink in
    % the final time step. The optional argument yDistRange
    % specifies the plot limits on the y axis in this plot, as a
    % multiplicative factor applied to the (automatically
    % determined) maximum. Default: yDistRange = 1e-20
    %
    % Clicking on a point in the heat-sink-magnitude plot produces
    % a new 1D plot with the parallel distribution at that time.
    %
    % Usage:
    %   HeatSink()
    %   HeatSink(yDistRange)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    yDistRange = 1e-20; %Range of source/distribution values (below 
                        %the maximum) to show.

    %%%%%%%%%%%%%%%%

    %Handle input
    switch nargin
        case 1
            %Nothing to do
        case 2
            if isnumeric(varargin{1}) && isscalar(varargin{1})
                yDistRange = varargin{1};
            else
                warning('The input must be a scalar.'); 
            end
        otherwise
            error('Invalid number of input arguments.');
    end

    oN = o.norse;
    figure(o.GetFigId('HeatSink'));
    clf;            

    %%% Heat sink magnitude evolution %%%
    subplot(2,1,1);

    %The signs are reversed here to get the correct appearance. The
    %time-advancement scheme changes it for the actual operator in
    %NORSE. The brackets are necessary since energyChangeMagnitudes
    %is a struct array.
    hs = zeros(1,8);
    mags = o.norse.energyChangeMagnitudes;
    nrm = mags.normalization;
    hs(1) = plot(oN.times,-[mags.sink]*nrm,'k','LineWidth',1.5,...
         'DisplayName','Heat sink/source','ButtonDownFcn',@o.PlotAtTimeX);
    hold on;
    hs(2) = plot(oN.times,-[mags.E]*nrm,'--',...
         'DisplayName','E-field acc.','ButtonDownFcn',@o.PlotAtTimeX);
    hs(3) = plot(oN.times,-[mags.C]*nrm,'--',...
         'DisplayName','Collisional redistribution','ButtonDownFcn',@o.PlotAtTimeX);
    hs(4) = plot(oN.times,-[mags.synch]*nrm,':',...
         'DisplayName','Synchrotron radiation reaction','ButtonDownFcn',@o.PlotAtTimeX);
    hs(5) = plot(oN.times,-[mags.tempChange]*nrm,'-.',...
         'DisplayName','Temperature changes','ButtonDownFcn',@o.PlotAtTimeX);             
    hs(6) = plot(oN.times,-[mags.densityChange]*nrm,'--c',...
         'DisplayName','Density changes','ButtonDownFcn',@o.PlotAtTimeX); 
    hs(7) = plot(oN.times,-[mags.correction]*nrm,'-.',...
         'DisplayName','Correction','ButtonDownFcn',@o.PlotAtTimeX); 
    hs(8) = plot(oN.times,-[mags.restriction]*nrm,':',...
         'DisplayName','Heat-loss-rate restriction','ButtonDownFcn',@o.PlotAtTimeX);

    hL = legend(hs,'Location','best');
    hL.Interpreter = 'latex';
    xlabel('$\tau$','Interpreter','latex');
    ylabel('$W/m^3$','Interpreter','latex');
    title('\bf{Energy addition/removal rate}','Interpreter','latex');



    %%% Final sink shape %%%
    subplot(2,1,2);
    if isempty(oN.heatSink)
        return
    end
    %Prepare the distribution and the sink
    allP = [-flipud(o.grid.p);o.grid.p];                   
    fOnGrid = o.grid.MapBigVectorToGrid(oN.f(:,end));
    f = [flipud(fOnGrid(:,1));fOnGrid(:,end)];

    sinkBig = (oN.dt*oN.heatSink.sink)*oN.f(:,end); %Apply the sink to f
    sinkOnGrid = o.grid.MapBigVectorToGrid(sinkBig);            
    sink = [flipud(sinkOnGrid(:,1));sinkOnGrid(:,end)];            

    %Determine plot limits
    yMax = 1.1*max([max(f),max(abs(sink))]);
    yMin = yDistRange*yMax;            
    xMaxId = min([find(fOnGrid(:,end)<yMin,1),o.grid.nP]);
    xMinId = min([find(fOnGrid(:,1)<yMin,1),o.grid.nP]);
    xMax = 1.1*o.grid.p(xMaxId);
    xMin = -1.1*o.grid.p(xMinId);

    %Plot it
    h1 = semilogy(allP,f,'-k','DisplayName','$F$');
    hold on;                        
    h2 = semilogy(allP,sink,'-r',...
                'DisplayName','$|\mathrm{d}\tau S_{\mathrm{h}}|$');
    semilogy(allP,-sink,':r');
    xlim([xMin,xMax])
    ylim([yMin,yMax])

    %Label the plot
    hL = legend([h1,h2],'Location','Best');
    hL.Interpreter = 'latex';
    xlabel('$p_{||}$','Interpreter','latex')            
    title('\bf{Heat sink/source, final time step}','Interpreter','latex')
end          
