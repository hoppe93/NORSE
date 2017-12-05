function out = Visualize(o,varargin)
    % Plots the spline and marks the specified data points. The
    % optional argument specifies the figure in which to plot. If
    % the parameter is constant in time, information is instead
    % printed to the console.
    %
    % Usage:
    %   Visualize()
    %   Visualize(hFig)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nPoints = 300; %to use for plotting the function value

    if o.isScalar
         fprintf('The parameter is constant, with the value: %.5g.\n',o.parameterSpline);
         out = []; %This is needed for handling the modified indexing
    else
        if nargin == 2
            hFig = figure(varargin{1});
        else
            hFig = figure(555);
        end                
        clf;

        dataPoints = o.parameterSpline.breaks; %Specified points
        pts = linspace(o.tMin,o.tMax,nPoints); 

        plot(pts,o.Evaluate(pts),'-k');
        hold on;
        plot(dataPoints,o.Evaluate(dataPoints),'*r');

        xlabel('Time');
        ylabel('Parameter value');
        xlim([o.tMin,o.tMax]);
        out = []; %This is needed for handling the modified indexing
    end            
end        
