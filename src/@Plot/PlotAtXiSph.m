function PlotAtXiSph(o,idT,~,eventdata)
    % Make a 1D plot at a given xi, specified by a click in a 2D 
    % plot in (p,xi) space.
    %
    % Usage: 
    %   Specify as ButtonDownFcn when plotting a quantity vs. time.
    %   The first argument must be passed using an anonymous
    %   funtion. Example:
    %       idT = 3;
    %       btnFcn = @(obj,event) PlotAtXiSph(idT,obj,event); 
    %       contourf(xi,p,f,'ButtonDownFcn',btnFcn);
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xi = eventdata.IntersectionPoint(1);
    [~,id] = min(abs(o.grid.xi-xi)); %Get the grid point closest to the clicked xi
    o.Dist1DAtXi(idT,id);
end 
