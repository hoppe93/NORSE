function PlotAtXiCyl(o,idT,~,eventdata)
    % Make a 1D plot at a given xi, specified by a click in a 2D 
    % plot in (p_para, p_perp) space.
    %
    % Usage: 
    %   Specify as ButtonDownFcn when plotting a quantity vs. time.
    %   The first argument must be passed using an anonymous
    %   funtion. Example:
    %       idT = 3;
    %       btnFcn = @(obj,event) PlotAtXiCyl(idT,obj,event); 
    %       contourf(pPara,pPerp,f,'ButtonDownFcn',btnFcn);
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    par = eventdata.IntersectionPoint(1);
    perp = eventdata.IntersectionPoint(2);
    xi = par/sqrt(par^2+perp^2);
    [~,id] = min(abs(o.grid.xi-xi)); %Get the grid point closest to the clicked xi
    o.Dist1DAtXi(idT,id);
end
