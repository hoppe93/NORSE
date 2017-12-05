function FunctionOfTime(o,vals,labelT,varargin)
    % Plots a vector as a function of time, with some convenient
    % labels and additions. The length of the vector must coincide
    % with the number of saved time step. 
    %
    % If an optional third argument is passed (and is true), the
    % x-axis label (\tau) is omitted. If a fourth argument is
    % passed, it determines the line style:
    %   0 = solid line 
    %   1 = dotted line
    %   2 = mark the data points (line style: '-*')            
    % 
    % Usage:
    %   FunctionOfTime(vals,yLabelText)
    %   FunctionOfTime(vals,yLabelText,suppressXLabel)
    %   FunctionOfTime(vals,yLabelText,suppressXLabel,lineStyle)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = o.norse.times;
    if nargin < 5 || ~varargin{2}
        plot(t,vals,'-b','ButtonDownFcn',@o.PlotAtTimeX);
    elseif varargin{2} == 1
        plot(t,vals,':b','ButtonDownFcn',@o.PlotAtTimeX);
    elseif varargin{2} == 2
        plot(t,vals,'-*b','ButtonDownFcn',@o.PlotAtTimeX);
    else
        error('Invalid line style argument.');
    end
    if ~(nargin >=4 && varargin{1})
        xlabel('$\tau$','Interpreter','latex');                
    end
    ylabel(labelT,'Interpreter','latex');
end
