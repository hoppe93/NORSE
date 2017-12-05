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
            o.collisionOperator = CLinear(o);
            if o.collisionOperatorMode == 1 %Full linear operator
                str = 'linear';
            else %Only test-particle operator
                str = 'linear, using only the test-particle term';                
            end
            switch o.timeAdvanceMode
                case 0
                    o.timeAdvance = TALinearDirect(o);
                case 1
                    error('This time-advance mode is not available yet for the linear collision operator.');
%                     o.timeAdvance = TAIterative(o);
                case 2
                    error('This time-advance mode is not available yet for the linear collision operator.');
%                     o.timeAdvance = TAAdaptive(o);
                case 3
                    error('This time-advance mode is not available yet for the linear collision operator.');
%                     o.timeAdvance = TAInductiveE(o);
                otherwise
                    error('Invalid time-advance mode');
            end            
        otherwise
            error('Invalid collision operator.')
    end  
    
    o.Print('   The collision operator is %s.\n',str);
end
