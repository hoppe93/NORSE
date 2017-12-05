function map = ColorMap(varargin)
    % Defines a color map that is suitable for both color and b&w
    % printing. 
    %
    % Usage:
    %   colormap(ColorMap())       -- Sets the colormap with 9 colors
    %   colormap(ColorMap(nColors))-- Sets the colormap with a
    %                                 specified number of colors
    %   colors = ColorMap(nColors) -- Returns an (nColors x 3) array
    %                                 of RGB triplets from the map
    %
    % The map was created by Gergely Papp. 
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    map = [ 0   0   0;...
           .15 .15 .5;...
           .3  .15 .75;...
           .6  .2  .50;...
            1  .25 .15;...
           .9  .5   0;...
           .9  .75 .1;...
           .9  .9  .5;...
            1   1   1];
    nColsInMap = size(map,1);

    if nargin>=1 
        nColors = varargin{1};
        if isscalar(varargin{1}) && isnumeric(varargin{1})
            %Interpolate the map to the desired number of colors
            colors = linspace(1,nColsInMap,nColors);
            map = interp1(map,colors);  
        else
            error('The argument must be a scalar');
        end
    end
end
