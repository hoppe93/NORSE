function AssembleSynchrotronOperator(o)
    % Builds the operator describing synchrotron-radiation-reaction
    % losses.
    %
    % Usage:
    %   AssembleSynchrotronOperator()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Build the operator
    if o.B.isScalar && o.B(0)==0
        %No synchrotron force -- set an empty operator
        o.synchrotronOperator = sparse(o.grid.matSize,o.grid.matSize);                
    else
        p = o.grid.pBig;
        xi = o.grid.xiBig;
        gamma = o.grid.gammaBig;
        oMinXi2 = 1-xi.*xi;

        %Helper function
        dia = @(v) spdiags(v,0,o.grid.matSize,o.grid.matSize);                

        dfdpTerm = dia(-oMinXi2.*gamma.*p)*o.grid.ddpMat;
        dfdxiTerm = dia(oMinXi2.*xi./gamma)*o.grid.ddxiMat;
        fTerm = dia( -oMinXi2./gamma.*(4*p.*p+2) - 2*xi.*xi./gamma );
        o.synchrotronOperator = -(dfdpTerm + dfdxiTerm + fTerm); 
                    %This is -R_S in the notation of the
                    %documentation. We change the sign as dictated
                    %by the time-advancement scheme
    end

    %Print some info
    if o.B.isScalar 
        if o.B(0) == 0 
            o.Print('   No synchrotron force is included.\n');
        elseif o.sigma.isScalar
            o.Print('   The magnetic field is %.2g T and sigma = %.3f.\n',...
                        o.B(0),o.sigma(0));
        else
            o.Print('   The magnetic field is %.2g T and sigma is time-dependent.\n',o.B(0));
        end                    
    else
        o.Print('   A time-dependent magnetic field is included.\n');
    end            
end 
