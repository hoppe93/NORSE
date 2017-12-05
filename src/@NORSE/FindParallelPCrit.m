function varargout = FindParallelPCrit(o,varargin)
    % Finds where dp/dt at xi=1 vanishes -- this is the end of the
    % separatrix. Optionally returns the sum of forces at all p on
    % the parallel grid (at xi=1) as the array dpdt. 
    % 
    % If the string 'inst' is passed as the first argument, the
    % current state of the distribution is used (specifically the
    % potential Pi from the Potentials object), and the second
    % argument is then the time. If the first argument is numeric,
    % it is interpreted as an index in the save array Pi. An
    % optional second argument determines whether to calculate the
    % force balance at xi=-1 or xi=1 (default).
    %
    % Usage: 
    %           pc = o.FindParallelPCrit('inst',t)
    %           pc = o.FindParallelPCrit(iteration)
    %    [pc,dpdt] = o.FindParallelPCrit(iteration)
    %    [pc,dpdt] = o.FindParallelPCrit(iteration,isAntiParallel)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    idXi = o.grid.nXi; %xi=1            
    switch nargin %o is the first argument, and will add +1 to the expected nargin values
        case 2
            iteration = varargin{1};
            EEc = o.EHat(o.times(iteration));
            Pi = o.Pi(:,iteration);                     
        case 3
            if strcmpi(varargin{1},'inst')
                Pi = o.potentials.Pi;
                EEc = o.EHat(varargin{2});
            else
                iteration = varargin{1};
                EEc = o.EHat(o.times(iteration));
                Pi = o.Pi(:,iteration); 
                if varargin{2}
                    idXi = 1; %xi=-1
                    EEc = -EEc; %E field now points towards bulk 
                end
            end    
        otherwise
            error('Wrong number of input arguments');
    end

    if nargout == 1 && EEc <= 1 
        %Don't do the calculation if the field is less than
        %critical in a standard call inside NORSE
        pc = Inf;
        dpdtOnGrid = zeros(size(o.grid.p));
    else
        %%%Find the critical momentum

        t0 = o.referenceTime;

        %%% Get the potential Pi, needed for the friction force
        Pi = o.grid.MapBigVectorToGrid(Pi);             
        PiAt1 = pchip(o.grid.p,Pi(:,idXi)); %Pi(xi=1) or Pi(xi=-1)
        dPidp = fnder(PiAt1);

        %%% Form the sum of forces on the parallel axis
        TNorm = o.lnLambdaBar(t0)/(o.Theta(t0)*o.kappa(t0)); 
        g2 = @(p) 1+p.*p;
        g = @(p) sqrt(g2(p));                                              
        dpdt = @(p) EEc + TNorm.*g(p) .* ppval(dPidp,p); %Synchrotron force vanishes at xi=1
        dpdtOnGrid = dpdt(o.grid.p);
        neg = dpdtOnGrid<0;

        if EEc<=1
            %No runaway region 
            pc = Inf;
        else
            %There is a runaway region -- find p_crit
            if ~any(neg)
            %Sum of forces positive everywhere -- everything runs away
            pc = 0;
        elseif ~any(dpdtOnGrid>0)
            %No runaway region on the grid
            pc = Inf;
        elseif sum(neg) < 2
            %Only a single point where force is negative. The next
            %point is the critical momentum
            idPc = find(neg,1)+1;
            pc = o.grid.p(idPc);
        else
            %General case -- find where force becomes positive
            guess = o.runawayRegion.pcs(end);
            if isinf(guess)
                guess = 1/sqrt(EEc-1);
            end                     
                pc = fzero(dpdt,guess);
                pc = min([pc,o.grid.p(end)]);
            end
        end
    end
    switch nargout
        case 0
            varargout = cell(1,1);
        case 1
            varargout = {pc};
        case 2
            varargout = {pc,dpdtOnGrid};
        otherwise
            error('Wrong number of output arguments');
    end
end
