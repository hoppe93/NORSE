function InitializeTAandC(o)
    % Initializes the proper TimeAdvance and CollisionOperator objects,
    % depending on the settings.
    %
    % Usage:
    %   InitializeTAandC()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    switch o.collisionOperatorMode
        case 0 %Nonlinear collision operator
            o.collisionOperator = CNonLinear(o);            
            str = 'nonlinear';
            switch o.timeAdvanceMode
                case 0
                    o.timeAdvance = TADirect(o);
                case 1
                    o.timeAdvance = TAIterative(o);
                case 2
                    o.timeAdvance = TAAdaptive(o);
                case 3
                    o.timeAdvance = TAInductiveE(o);
                otherwise
                    error('Invalid time-advance mode');
            end
        case {1,2} %Linear collision operator
			error('The linear collision operator is not available yet.');                       
        otherwise
            error('Invalid collision operator.')
    end  
    
    o.Print('   The collision operator is %s.\n',str);
end
