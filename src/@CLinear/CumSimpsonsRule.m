function cumInt = CumSimpsonsRule(f,p)     
    % Calculate the cumulative integral of a function f(p) over the grid p.
    % The integral on each interval in p is calculated using Simpson's rule
    % with 4 subintervals. f must be a function handle taking one
    % vectorized argument.
    %
    % Usage:
    %   cumInt = CumSimpsonsRule(f,p)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Define four shifted grids
    dp = diff(p);    
    p1 = p(1:end-1);
    p2 = p1+0.25*dp;
    p3 = p1+0.5*dp;
    p4 = p1+0.75*dp;    
    
    %Evaluate the function on the grids and put the results in an
    %(ny x 5) matrix including the quadrature weights 
    fp = f(p);
    intMat=[fp(1:end-1),...
             4*f(p2),...
             2*f(p3),...
             4*f(p4),...                  
             fp(2:end)];
         
    %Perform the integeration as a summation over rows, then calculate the
    %cumulative integral
    cumInt = 1/12*dp.*sum(intMat,2);
    cumInt = [0;cumsum(cumInt)];   
end
