function PlotAtTimeY(o,~,eventdata)
    % Plots the 1D distribution at the time selected with the
    % mouse, when the time coordinate is on the y axis.
    %
    % Usage: 
    %   Specify as ButtonDownFcn when plotting a quantity vs. time.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = eventdata.IntersectionPoint(2);            
    [~,id] = min(abs(o.norse.times-t)); %Get the saved point closest to the clicked time
    o.Dist1D(id);
end 
