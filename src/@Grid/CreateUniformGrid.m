function [x,w,ddx,d2dx2] = CreateUniformGrid(o,n,xMin,xMax)
    % Generates a uniform grid, together with quadrature weights
    % and the first and second finite-diffrences. The quadrature
    % weights can be chosen according to either the trapezoid rule
    % or a composite Simpson's rule (which should be more
    % accurate). Stencils with 3,5 or 7 points (giving 2nd, 4th or
    % 6th order finite differences) are available. Central
    % differences are used, except near the endpoints where a
    % gradual transition to one-sided differences is used to
    % maintain the discretization order.
    %
    % Usage:
    %   [x,weights,ddx,d2dx2] = CreateUniformGrid(nPoints,xMin,xMax)
    %
    % Dimensions:
    %   x, weights: [nPoints,1]            
    %   ddx, d2dx2: [nPoints,nPoints]
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Validate input:
    if ~isfinite(n)
        error('n of points must be finite')
    end
    if ~isfinite(xMin)
        error('xMin must be finite')
    end
    if ~isfinite(xMax)
        error('xMax must be finite')
    end
    if numel(n) ~= 1
        error('n must be a single integer')
    end
    if numel(xMin) ~= 1
        error('xMin must be a single number')
    end
    if numel(xMax) ~= 1
        error('xMax must be a single number')
    end    
    if xMax < xMin
        error('xMax must be larger than xMin.')
    end
    if xMax == xMin
        error('xMax cannot equal xMin.')
    end
    if n<o.stencil
        error('n must be at least the number of points of the stencil.')
    end
    if n ~= round(n)
        error('n must be an integer.')
    end

    % Create grid
    x   = linspace(xMin,xMax,n)';
    dx  = x(2)-x(1);
    dx2 = dx*dx;

    % Create integration weights
    switch lower(o.quadrature)
        case 'trapezoid'
            w      = ones(size(x));    
            w(1)   = 0.5;
            w(end) = 0.5;    
            w      = w*dx;
        case 'simpson'
            %Use a composite Simpson's rule. Use simpler rules in
            %the (unlikely) case that we have fewer than 8 points.
            switch n            
                case 3
                    w = [1,4,1]'/6; %Simpson's rule
                case 4
                    w = 3/8*[1,3,3,1]'; %Simpson's 3/8 rule
                case 5
                    w = [1,4,2,4,1]'/3; %Composite Simpson's rule
                case 6
                    w = 5/288*[19,75,50,50,75,19]'; %6-point Newton-Coates rule
                case 7
                    w = [1,4,2,4,2,4,1]'/3; %Composite Simpson's rule
                otherwise
                    w = 1/48*[17,59,43,49,48*ones(1,(n-8)),49,43,59,17]';
            end            
            w = w*dx;            
        otherwise
            error('Invalid quadrature.');
    end

    % Create differentiation matrices
    switch o.stencil
        case 3  
            %%% ddx
            % Set the interior points of the differentiation matrix
            % using a central finite difference:                    
            ddx = sparse( (  diag(ones(n-1,1),1) ...
                           - diag(ones(n-1,1),-1) )...
                         /(2*dx) );

            %Handle endpoints (one-sided differences):
            ddx(1,1) = -1.5/dx;
            ddx(1,2) = 2/dx;
            ddx(1,3) = -0.5/dx;

            ddx(end,end)   = 1.5/dx;
            ddx(end,end-1) = -2/dx;
            ddx(end,end-2) = 0.5/dx;


            %%% d2dx2
            %Interior points
            d2dx2 = sparse( ( - 2*diag(ones(n,1),0) ...
                              + diag(ones(n-1,1),1)...
                              + diag(ones(n-1,1),-1) )...
                            /dx2 );

            %Endpoints
            d2dx2(1,1) = 1/dx2;
            d2dx2(1,2) = -2/dx2;
            d2dx2(1,3) = 1/dx2;

            d2dx2(end,end)   = 1/dx2;
            d2dx2(end,end-1) = -2/dx2;
            d2dx2(end,end-2) = 1/dx2;            
        case 5  
            %%% ddx
            %Interior points
            v1  = ones(n,1);
            ddx = (  - spdiags((1/6)*v1,2,n,n)...
                     + spdiags((4/3)*v1,1,n,n)...
                     - spdiags((4/3)*v1,-1,n,n)... 
                     + spdiags((1/6)*v1,-2,n,n) )...
                   /(2*dx);

            %Endpoints
            ddx(1,1) = -25/(12*dx);
            ddx(1,2) = 4/(dx);
            ddx(1,3) = -3/dx;
            ddx(1,4) = 4/(3*dx);
            ddx(1,5) = -1/(4*dx);

            ddx(2,1) = -1/(4*dx);
            ddx(2,2) = -5/(6*dx);
            ddx(2,3) = 3/(2*dx);
            ddx(2,4) = -1/(2*dx);
            ddx(2,5) = 1/(12*dx);

            ddx(end,end)   = 25/(12*dx);
            ddx(end,end-1) = -4/(dx);
            ddx(end,end-2) = 3/dx;
            ddx(end,end-3) = -4/(3*dx);
            ddx(end,end-4) = 1/(4*dx);

            ddx(end-1,end)   = 1/(4*dx);
            ddx(end-1,end-1) = 5/(6*dx);
            ddx(end-1,end-2) = -3/(2*dx);
            ddx(end-1,end-3) = 1/(2*dx);
            ddx(end-1,end-4) = -1/(12*dx);


            %%% d2dx2
            %Interior points
            d2dx2 = (   - spdiags((1/12)*v1,2,n,n)...
                        + spdiags((4/3)*v1,1,n,n)...
                        + spdiags((-5/2)*v1,0,n,n)...
                        + spdiags((4/3)*v1,-1,n,n) ...
                        - spdiags((1/12)*v1,-2,n,n) )...
                      /(dx2);

            %Endpoints
            d2dx2(1,1) = 35/(12*dx2);
            d2dx2(1,2) = -26/(3*dx2);
            d2dx2(1,3) = 19/(2*dx2);
            d2dx2(1,4) = -14/(3*dx2);
            d2dx2(1,5) = 11/(12*dx2);

            d2dx2(2,1) = 11/(12*dx2);
            d2dx2(2,2) = -5/(3*dx2);
            d2dx2(2,3) = 1/(2*dx2);
            d2dx2(2,4) = 1/(3*dx2);
            d2dx2(2,5) = -1/(12*dx2);

            d2dx2(end,end)   = 35/(12*dx2);
            d2dx2(end,end-1) = -26/(3*dx2);
            d2dx2(end,end-2) = 19/(2*dx2);
            d2dx2(end,end-3) = -14/(3*dx2);
            d2dx2(end,end-4) = 11/(12*dx2);

            d2dx2(end-1,end-0) = 11/(12*dx2);
            d2dx2(end-1,end-1) = -5/(3*dx2);
            d2dx2(end-1,end-2) = 1/(2*dx2);
            d2dx2(end-1,end-3) = 1/(3*dx2);
            d2dx2(end-1,end-4) = -1/(12*dx2);            
        case 7  
            v1 = ones(n,1);

            %%% ddx
            %Interior points
            ddx = (    spdiags((1/60)*v1,3,n,n)...
                     + spdiags((-3/20)*v1,2,n,n) ...
                     + spdiags((3/4)*v1,1,n,n)...
                     + spdiags((-3/4)*v1,-1,n,n)...
                     + spdiags((3/20)*v1,-2,n,n)...
                     + spdiags((-1/60)*v1,-3,n,n) )...
                   /dx;


            %Endpoints
            ddx(1,1) = -49/(20*dx);
            ddx(1,2) = 6/dx;
            ddx(1,3) = -15/(2*dx);
            ddx(1,4) = 20/(3*dx);
            ddx(1,5) = -15/(4*dx);
            ddx(1,6) = 6/(5*dx);
            ddx(1,7) = -1/(6*dx);

            ddx(2,1) = -1/(6*dx);
            ddx(2,2) = -77/(60*dx);
            ddx(2,3) = 5/(2*dx);
            ddx(2,4) = -5/(3*dx);
            ddx(2,5) = 5/(6*dx);
            ddx(2,6) = -1/(4*dx);
            ddx(2,7) = 1/(30*dx);

            ddx(3,1) = 1/(30*dx);
            ddx(3,2) = -2/(5*dx);
            ddx(3,3) = -7/(12*dx);
            ddx(3,4) = 4/(3*dx);
            ddx(3,5) = -1/(2*dx);
            ddx(3,6) = 2/(15*dx);
            ddx(3,7) = -1/(60*dx);



            ddx(end,end)   = 49/(20*dx);
            ddx(end,end-1) = -6/dx;
            ddx(end,end-2) = 15/(2*dx);
            ddx(end,end-3) = -20/(3*dx);
            ddx(end,end-4) = 15/(4*dx);
            ddx(end,end-5) = -6/(5*dx);
            ddx(end,end-6) = 1/(6*dx);

            ddx(end-1,end)   = 1/(6*dx);
            ddx(end-1,end-1) = 77/(60*dx);
            ddx(end-1,end-2) = -5/(2*dx);
            ddx(end-1,end-3) = 5/(3*dx);
            ddx(end-1,end-4) = -5/(6*dx);
            ddx(end-1,end-5) = 1/(4*dx);
            ddx(end-1,end-6) = -1/(30*dx);

            ddx(end-2,end)   = -1/(30*dx);
            ddx(end-2,end-1) = 2/(5*dx);
            ddx(end-2,end-2) = 7/(12*dx);
            ddx(end-2,end-3) = -4/(3*dx);
            ddx(end-2,end-4) = 1/(2*dx);
            ddx(end-2,end-5) = -2/(15*dx);
            ddx(end-2,end-6) = 1/(60*dx);





            %%% d2dx2
            %Interior points
            d2dx2 = (     spdiags((1/90)*v1,3,n,n)...
                        + spdiags((-3/20)*v1,2,n,n)...
                        + spdiags((3/2)*v1,1,n,n)...
                        + spdiags((-49/18)*v1,0,n,n)...
                        + spdiags((3/2)*v1,-1,n,n)...
                        + spdiags((-3/20)*v1,-2,n,n)...
                        + spdiags((1/90)*v1,-3,n,n)  )...
                     /dx2;


            %Endpoints     
            d2dx2(1,1) = 203/(45*dx2);
            d2dx2(1,2) = -87/(5*dx2);
            d2dx2(1,3) = 117/(4*dx2);
            d2dx2(1,4) = -254/(9*dx2);
            d2dx2(1,5) = 33/(2*dx2);
            d2dx2(1,6) = -27/(5*dx2);
            d2dx2(1,7) = 137/(180*dx2);

            d2dx2(2,1) = 137/(180*dx2);
            d2dx2(2,2) = -49/(60*dx2);
            d2dx2(2,3) = -17/(12*dx2);
            d2dx2(2,4) = 47/(18*dx2);
            d2dx2(2,5) = -19/(12*dx2);
            d2dx2(2,6) = 31/(60*dx2);
            d2dx2(2,7) = -13/(180*dx2);

            d2dx2(3,1) = -13/(180*dx2);
            d2dx2(3,2) = 19/(15*dx2);
            d2dx2(3,3) = -7/(3*dx2);
            d2dx2(3,4) = 10/(9*dx2);
            d2dx2(3,5) = 1/(12*dx2);
            d2dx2(3,6) = -1/(15*dx2);
            d2dx2(3,7) = 1/(90*dx2);



            d2dx2(end,end)   = 203/(45*dx2);
            d2dx2(end,end-1) = -87/(5*dx2);
            d2dx2(end,end-2) = 117/(4*dx2);
            d2dx2(end,end-3) = -254/(9*dx2);
            d2dx2(end,end-4) = 33/(2*dx2);
            d2dx2(end,end-5) = -27/(5*dx2);
            d2dx2(end,end-6) = 137/(180*dx2);

            d2dx2(end-1,end)   = 137/(180*dx2);
            d2dx2(end-1,end-1) = -49/(60*dx2);
            d2dx2(end-1,end-2) = -17/(12*dx2);
            d2dx2(end-1,end-3) = 47/(18*dx2);
            d2dx2(end-1,end-4) = -19/(12*dx2);
            d2dx2(end-1,end-5) = 31/(60*dx2);
            d2dx2(end-1,end-6) = -13/(180*dx2);

            d2dx2(end-2,end)   = -13/(180*dx2);
            d2dx2(end-2,end-1) = 19/(15*dx2);
            d2dx2(end-2,end-2) = -7/(3*dx2);
            d2dx2(end-2,end-3) = 10/(9*dx2);
            d2dx2(end-2,end-4) = 1/(12*dx2);
            d2dx2(end-2,end-5) = -1/(15*dx2);
            d2dx2(end-2,end-6) = 1/(90*dx2);            
        otherwise
            error('Invalid stencil.');
    end    
end 
