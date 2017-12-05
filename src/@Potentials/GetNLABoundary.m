function NLA = GetNLABoundary(o,l,a) 
    % Generate a vector NLA with the p' dependence of N_l,a at pMax.
    % This can then be used to integrate over p'. 
    %
    % Usage:
    %   NLA = GetNLABoundary(l,a)
    %
    % l is the Legendre-mode index and a is the parameter in the
    % differential operator L_a.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    switch a
        case '0'
            yL0 = o.jAndYMatrices.yL0(:,l+1);
            jL0 = o.jAndYMatrices.jL0(:,l+1);
            NLA = yL0(end)*jL0'; %This gives the last row of the 
                                 %full NLA matrix yL0*jL0'
        case '02'
            yL0 = o.jAndYMatrices.yL0(:,l+1);
            jL02 = o.jAndYMatrices.jL02(:,l+1);
            yL02 = o.jAndYMatrices.yL02(:,l+1);
            jL2 = o.jAndYMatrices.jL2(:,l+1);
            NLA = yL0(end)*jL02' + yL02(end)*jL2';
        case '022'
            yL0 = o.jAndYMatrices.yL0(:,l+1);
            jL022 = o.jAndYMatrices.jL022(:,l+1);
            yL02 = o.jAndYMatrices.yL02(:,l+1);
            jL22 = o.jAndYMatrices.jL22(:,l+1);
            yL022 = o.jAndYMatrices.yL022(:,l+1);
            jL2 = o.jAndYMatrices.jL2(:,l+1);
            NLA = yL0(end)*jL022' + yL02(end)*jL22' + yL022(end)*jL2';
        case '1'
            yL1 = o.jAndYMatrices.yL1(:,l+1);
            jL1 = o.jAndYMatrices.jL1(:,l+1);
            NLA = yL1(end)*jL1';
        case '11'
            yL1 = o.jAndYMatrices.yL1(:,l+1);
            jL11 = o.jAndYMatrices.jL11(:,l+1);
            yL11 = o.jAndYMatrices.yL11(:,l+1);
            jL1 = o.jAndYMatrices.jL1(:,l+1);
            NLA = yL1(end)*jL11' + yL11(end)*jL1';
        otherwise
            error('Invalid combination of indices a.');
    end
    NLA(isinf(NLA)) = realmax;
    NLA(isnan(NLA)) = 0;
end        
