function Initialize1DTimeEvolutionPlot(o)
    % Initializes a plot for showing a parallel cut of the
    % distribution during runtime.
    %
    % Usage:
    %   Initialize1DTimeEvolutionPlot()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    lw = 2; %Linewidth in the plot

    %%%%%%%%

    warning off MATLAB:Axes:NegativeDataInLogAxis                

    figure(o.figOffset + 1);
    clf;

    F = o.grid.MapBigVectorToGrid(o.f(:,1));
    F = [flipud(F(:,1)); F(:,end)];
    Fpos = 1e-100*ones(size(F));
    Fneg = Fpos;
    Fpos(F>0) = F(F>0);
    Fneg(F<=0) = F(F<=0);            

    p = [-flipud(o.grid.p); o.grid.p];
    axis([min(p)/2 max(p) 1e-25 1])
    semilogy(p,F,'--k','linewidth',lw)            
    hold on
    o.timeEvoLine{1} = semilogy(p,Fpos,'-k','linewidth',lw);
    hold on
    if any(F<=0)
        o.timeEvoLine{2} = semilogy(p,-Fneg,':k','linewidth',lw);
    else
        o.timeEvoLine{2} = semilogy(p,zeros(size(Fneg)),':k','linewidth',lw);
    end

    xlabel('$p_{||}$','fontsize',14,'Interpreter','latex');
    ylabel('$F$ (at $\xi\!=\!1$)','fontsize',14,'Interpreter','latex');
    leg = legend('Initial','$F_{\mathrm{pos}}(t)$','$F_{\mathrm{neg}}(t)$');
    set(gca,'fontsize',12,'linewidth',lw,'ytick',10.^(-25:5:0));
    set(leg,'linewidth',1.5,'Interpreter','latex');
end
