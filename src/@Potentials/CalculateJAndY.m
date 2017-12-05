function CalculateJAndY(o)
    % Calculate arrays containing j_l[k]* and y_l[k]* for all grid
    % points and all Legendre modes. The notation is that of Braams
    % & Karney [PoF B, 1 (1989), 1355] and the implementation is
    % based on the algorithm in their Appendix 7.
    %
    % The arrays are saved in the class property jAndYMatrices; a
    % struct with the fields:
    %   jL0, jL1, jL2
    %   yL0, yL1, yL2
    %   jL02, yL02
    %   jL11, yL11
    %   jL22, yL22
    %   jL022, yL022  
    %
    % Each field f=f(p,l) has the structure
    %   [f(p1,0), ..., f(p1,nL-1)]
    %   [ ...   , ...,   ...  ]
    %   [f(pN,0), ..., f(pN,nL-1)]
    %
    %   
    % To obtain the j_l[1]a for all l>=0, the routine uses forward
    % recursion for large p (where it works excellently), and
    % backward recursion close to p=0, where the forward method is
    % unreliable. The value of p for which the transition between
    % the two methods occurs is l-dependent, for optimal
    % performance. j_l[1]a for l<0 is calculated using backward
    % recursion for all p. The y_l[1]a, j_l[k]* and y_l[k]* are
    % then calculated from the series of j_l[1]a using simple
    % recurrence relations. 
    %
    % Usage:
    %   CalculateJAndY()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    delta = 1e-14; %Tolerance for the backward recursion
    splitDependence = @(l) 0 + l*0.9/25; 
                %Defines the dependence of the split point between
                %forward and backward recursion on l.

    %%%%%%%

    L = (o.norse.nL-1+2); %nL includes the l=0 mode. We also need two 
                          %extra L's for the higher order functions
    pFull = o.grid.p;
    sJY = struct(); %Struct for saving all the arrays

    % Let's calculate with forward recursion for all points (the
    % extra cost is small), but only with backward recursion for
    % points close to 0.
    splitPoint = splitDependence(L); 
    idSplit = find(pFull>splitPoint,1);
    if isempty(idSplit)
        idSplit = numel(pFull);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Do the backward recursion for the low-p part %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = pFull(1:idSplit);
    p2 = p.*p;
    gamma = sqrt(p.^2+1); 

    %%% Calculate L'
    pMax = p(end);
    LPrim = L + ceil(log(delta)/(2*log(pMax/(sqrt(pMax^2+1)+1))));        

    %%% Loop over l and determine j recursively %%%%%%%%%%%%%%%%%%%
    % First loop through the first L'-L iterations, to get to the
    % ones we actually want to save. In fact, let's go two steps
    % further so that the information in jLPlus1 and jLPlus2 can be
    % directly put into the array we want to save. j here is
    % actually \tilde{j}, which is more numerically stable.
    jLPlus2 = 2./(gamma+1); 
    jLPlus1 = ones(size(gamma));
    for l = LPrim:-1:(L-1)
        jL = gamma.*jLPlus1 - l*l/(4*l*l-1)*p2.*jLPlus2; 
                 %Eq.(4) with a=0. This gives overflow for large p!
        jLPlus2 = jLPlus1;
        jLPlus1 = jL;
    end

    % Now do the ones we need to save. We will populate the array
    % from the right.
    jL0(:,L+3) = jLPlus2;
    jL0(:,L+2) = jLPlus1;
    for l = (L-2):-1:-2
        iCol = l+3;
        jL0(:,iCol) = gamma.*jL0(:,iCol+1) ...
                      -  l*l/(4*l^2-1)*p2.*jL0(:,iCol+2);
    end

    %Renormalize from \tilde{j} to j
    iMinus1 = 2;
    jL0 = jL0./(jL0(:,iMinus1)*ones(1,(L+3)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Calculate jl,1
    jL1 = zeros(size(jL0));
    %Get initial conditions from jl,0 using Eq. 6
    jL1(:,L+3) = 1/(L+1) * ( (2*L+1)*jL0(:,end-1) - L*gamma.*jL0(:,end) ); %Eq 6 with l=L, a=1
    jL1(:,L+2) = 1/L * ( (2*L-1)*jL0(:,end-2) - (L-1)*gamma.*jL0(:,end-1) ); %Eq 6 with l=L-1, a=1
    for l = (L-2):-1:-2
        iCol = l+3;
        jL1(:,iCol) = gamma.*jL1(:,iCol+1) - (l-1)*(l+1)/(4*l^2-1)*p2.*jL1(:,iCol+2); %Eq.(4) with a=1 
    end

    %%% Calculate jl,2
    jL2 = zeros(size(jL0));
    %Get initial conditions from jl,1 using Eq. 6
    jL2(:,L+3) = 1/(L+2) * ( (2*L+1)*jL1(:,end-1) - (L-1)*gamma.*jL1(:,end) ); %Eq 6 with l=L, a=2
    jL2(:,L+2) = 1/(L+1) * ( (2*L-1)*jL1(:,end-2) - (L-2)*gamma.*jL1(:,end-1) ); %Eq 6 with l=L-1, a=2
    for l = (L-2):-1:-2
        iCol = l+3;
        jL2(:,iCol) = gamma.*jL2(:,iCol+1) - (l-2)*(l+2)/(4*l^2-1)*p2.*jL2(:,iCol+2); %Eq.(4) with a=2 
    end

    jL0Low = jL0;
    jL1Low = jL1;
    jL2Low = jL2;



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Do the forward recursion for the high-p part %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = pFull;
    p2 = p.*p;
    pInv2 = 1./p2;
    gamma = o.grid.gamma;
    sigma = asinh(p);    

    %%% Loop over l and determine j recursively %%%%%%%%%%%%%%%%%%%   

    %%% Calculate jL,0
    jL0 = zeros(numel(p),L+3);
    j010 = sigma./p; %Initial conditions from the known analytical 
                     %expressions. \tilde{j}=j for l=0.
    j110 = ((gamma.*sigma-p)./p2) * 3./p; %Normalized (\tilde{j})
    jL0(:,3) = j010;
    jL0(:,4) = j110;

    % We will populate the array from the middle to the right and
    % then from the middle to the left.
    % The l>0...
    for l = 2:L
        iCol = l+3;
        jL0(:,iCol) = (4*l*l-1)/(l*l)*pInv2 .* ...
                      ( gamma.*jL0(:,iCol-1)-jL0(:,iCol-2) ); %Eq 7 with a=0
    end    
    % ...and now the l<0 (backward recursion)
    for l = -1:-1:-2
        iCol = l+3;
        jL0(:,iCol) = gamma.*jL0(:,iCol+1) ...
                      - l*l/(4*l^2-1)*p2.*jL0(:,iCol+2); %Eq.(4) with a=0
    end


    %%% Calculate jl,1    
    jL1 = zeros(size(jL0));
    j011 = ones(size(p)); 
    j111 = (0.5*(p.*gamma-sigma)./p2) * 3./p;
    jL1(:,3) = j011;
    jL1(:,4) = j111;    
    for l = 2:L
        iCol = l+3;
        jL1(:,iCol) = (4*l*l-1)/((l-1)*(l+1))*pInv2 .* ...
                      ( gamma.*jL1(:,iCol-1)-jL1(:,iCol-2) ); %Eq 7 with a=1
    end        
    for l = -1:-1:-2
        iCol = l+3;
        jL1(:,iCol) = gamma.*jL1(:,iCol+1) - ...
                      (l-1)*(l+1)/(4*l*l-1)*p2.*jL1(:,iCol+2); %Eq.(4) with a=1
    end


    %%% Calculate jl,2    
    jL2 = zeros(size(jL0));
    j012 = gamma;
    j112 = ones(size(p));
    jL2(:,3) = j012;
    jL2(:,4) = j112;
    jL2(:,5) = 0.25 * ( 5*jL1(:,4) - gamma.*jL1(:,5) ); 
            %Eq 6 with l=2, a=2. We need this for a well behaved recursion.     
    for l = 3:L
        iCol = l+3;
        jL2(:,iCol) = (4*l*l-1)/((l-2)*(l+2))*pInv2 .* ...
                      ( gamma.*jL2(:,iCol-1)-jL2(:,iCol-2) ); %Eq 7 with a=2
    end        
    for l = -1:-1:-2
        iCol = l+3;
        jL2(:,iCol) = gamma.*jL2(:,iCol+1) - ...
                      (l-2)*(l+2)/(4*l*l-1)*p2.*jL2(:,iCol+2); %Eq.(4) with a=2
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Combine the low and high-p parts %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Substitue in results from backward calculation for points
    % close to p=0.
    for l = 0:L
        iCol = l+3;
        id = find(pFull>splitDependence(l),1);
        if isempty(id)
            id = numel(pFull);
        end
        jL0(1:id,iCol) = jL0Low(1:id,iCol); 
        jL1(1:id,iCol) = jL1Low(1:id,iCol); 
        jL2(1:id,iCol) = jL2Low(1:id,iCol); 
    end

    %%% Renormalize from \tilde{j} to {j}
    % We need to form the multiplication matrix p^l/(2*l+1)!! in a
    % loop to avoid overflow.
    reNoMat(size(jL0,1),size(jL0,2)) = 0; %reserve memory
    reNoMat(:,3) = ones(size(jL0,1),1); %The l=0 case is just 1    
    for l = 1:L %l>0
        iCol = l+3;
        reNoMat(:,iCol) = reNoMat(:,iCol-1).*p/(2*l+1);         
    end
    for l = -1:-1:-2 %l<0
        iCol = l+3;
        reNoMat(:,iCol) = (2*l+3)*reNoMat(:,iCol+1)./p;         
    end
    jL0Tilde = jL0;
    jL1Tilde = jL1;
    jL2Tilde = jL2;
    jL0 = reNoMat.*jL0;
    jL1 = reNoMat.*jL1;
    jL2 = reNoMat.*jL2;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Compute yLa %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Calculate jLa for l<0 using backward recursion %%%%%%%%%%%%            
    jL0Neg(numel(p),L+2) = 0;
    jL1Neg(numel(p),L+2) = 0;
    jL2Neg(numel(p),L+2) = 0;

    % Take values for -a<=l<=0 from previous calculation, and
    % define initial values for the backward recursion. Indexing is
    % iCol=-l+2.

    % l=1 -- we need this to compute the higher order modes
    oldId = 1+3;
    jL0Neg(:,1) = jL0Tilde(:,oldId);    
    jL1Neg(:,1) = jL1Tilde(:,oldId);    
    jL2Neg(:,1) = jL2Tilde(:,oldId);   

    % l=0
    oldId = 0+3;
    jL0Neg(:,2) = jL0Tilde(:,oldId);    
    jL1Neg(:,2) = jL1Tilde(:,oldId);    
    jL2Neg(:,2) = jL2Tilde(:,oldId);        

    % l=-1    
    oldId = -1+3;
    jL0Neg(:,3) = 1; %j_{-1,0}
    jL1Neg(:,3) = jL1Tilde(:,oldId); %j_{-1,1}
    jL2Neg(:,3) = jL2Tilde(:,oldId); %j_{-1,2}

    % l=-2
    oldId = -2+3;
    jL1Neg(:,4) = 1; %j_{-2,1}
    jL2Neg(:,4) = jL2Tilde(:,oldId); %j_{-2,2}    

    %l=-3
    jL2Neg(:,5) = 1; %j_{-3,2}


    % --- Do backward recursion ---
    % First, manually bring the a=0,1 cases to l=-3
    jL0Neg(:,4) = gamma.*jL0Neg(:,3); %(A29) with l-2=-2,a=0
    jL0Neg(:,5) = gamma.*jL0Neg(:,4) ...
                  - 1/3*p2.*jL0Neg(:,3); %(A29) with l-2=-3,a=0
    jL1Neg(:,5) = gamma.*jL1Neg(:,4); %(A29) with l-2=-3,a=1    

    % Recursively calculate the values for the remaining l
    for il = -4:-1:-L%-(L-2)
        l = -il+2;
        denom = (2*il+3)*(2*il+5);
        jL0Neg(:,l) = gamma.*jL0Neg(:,l-1) - ...
            (il+2)*(il+2)/denom*p2.*jL0Neg(:,l-2); %(A29) with l-2=il, a=0
        jL1Neg(:,l) = gamma.*jL1Neg(:,l-1) -...
            (il+1)*(il+3)/denom*p2.*jL1Neg(:,l-2); %(A29) with l-2=il, a=1
        jL2Neg(:,l) = gamma.*jL2Neg(:,l-1) - ...
            il*(il+4)/denom*p2.*jL2Neg(:,l-2); %(A29) with l-2=il, a=2
    end

    %Renormalize \tilde{j} to j
    reNoMat = zeros(size(jL0Neg)); %reserve memory
    reNoMat(:,1) = p/3; %The l=1 case
    reNoMat(:,2) = 1; %The l=0 case is just 1 
    for l = -1:-1:-L %l<0
        iCol = -l+2;
        reNoMat(:,iCol) = (2*l+3)*reNoMat(:,iCol-1)./p;         
    end
    jL0Neg = reNoMat.*jL0Neg;
    jL1Neg = reNoMat.*jL1Neg;
    jL2Neg = reNoMat.*jL2Neg;

    %%% Now calculate the yLa using Eq. (A8) %%%%%%%%%%%%%%%%%%%%%%
    yL0(numel(p),L+1) = 0;
    yL1(numel(p),L+1) = 0;
    yL2(numel(p),L+1) = 0;
    preFac = -1;
    for l = -2:L-1  
        iColJ = l+3; 
        yL0(:,iColJ) = preFac*jL0Neg(:,iColJ);
        yL1(:,iColJ) = preFac*jL1Neg(:,iColJ);
        yL2(:,iColJ) = preFac*jL2Neg(:,iColJ);
        preFac = -1*preFac;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Compute the higher order functions %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Calculate jl[k]* using Eq. (A26)
    %The jLa contain -2:L, but the jL* should contain only 0:L-2
    %(the extra two l's are only used for the formulas below)
    pHalf = 0.5*p*ones(1,L-1);
    p2Half = 0.5*p2*ones(1,L-1);
    sJY.jL02 = pHalf.*jL1(:,4:(end-1)); 
    sJY.jL11 = pHalf.*jL0(:,4:(end-1));
    sJY.jL22 = pHalf.*jL1(:,4:(end-1)) + p2Half.*jL0(:,5:end);
    sJY.jL022 = 0.25*p2Half.*jL0(:,5:end);


    %%% Calculate yl[k]* using Eq. (A26) (together with A8)
    %The yLa contain -2:(L-1), but the yL* should contain only 0:(L-1)
    pHalf = -pHalf;
    sJY.yL02 = pHalf.*yL1(:,2:(end-2));
    sJY.yL11 = pHalf.*yL0(:,2:(end-2));
    sJY.yL22 = pHalf.*yL1(:,2:(end-2)) + p2Half.*yL0(:,1:(end-3));
    sJY.yL022 = 0.25*p2Half.*yL0(:,1:(end-3));   


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %Trim the jLa and yLa arrays    
    sJY.jL0 = jL0(:,3:end-2);
    sJY.jL1 = jL1(:,3:end-2);
    sJY.jL2 = jL2(:,3:end-2);
    sJY.yL0 = yL0(:,3:end-1);
    sJY.yL1 = yL1(:,3:end-1);
    sJY.yL2 = yL2(:,3:end-1);

    o.jAndYMatrices = sJY;
end
