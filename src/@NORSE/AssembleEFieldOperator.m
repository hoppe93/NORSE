function AssembleEFieldOperator(o)
    % Builds the operator describing the electric field
    % acceleration.
    %
    % Usage: 
    %   AssembleEFieldOperator()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sz = o.grid.matSize;

    if (o.EHat.isScalar && o.EHat(1)==0) 
        %No electric field -- define an empty operator
        o.EFieldOperator = sparse(sz,sz);                
    else
        ddpTerm  = spdiags(o.grid.xiBig,0,sz,sz)*o.grid.ddpMat;
        ddxiTerm = spdiags((1-o.grid.xiBig.^2)./o.grid.pBig,0,sz,sz)*o.grid.ddxiMat;
        ddxiTerm(end,:)  = 0; %There is no xi-dependence at p=0 
        o.EFieldOperator = -ddpTerm - ddxiTerm; %Change the sign, to make the tail point to the right                
    end

    %Print some info about E and Z
    if o.EHat.isScalar 
        if o.EHat(0) == 0 
            o.Print('   No electric field is applied. Z_eff = %.1f.\n',o.Z(0));
        else
            o.Print('   The electric field is %.3g E/E_D (%.3g E/E_c) and Z_eff = %.1f.\n',...
                        o.Theta(0)*o.EHat(0),o.EHat(0),o.Z(0));
        end                    
    else
        o.Print('   A time-dependent electric field is applied. Z_eff = %.1f.\n',o.Z(0));
    end            
end 
