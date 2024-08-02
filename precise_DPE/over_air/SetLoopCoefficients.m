function LoopFilter = SetLoopCoefficients( bandwidth, order, Tc, initialValue )

LoopFilter = struct('inputValues', zeros(1, 3), ...
                    'outputValues', zeros(1, 2), ...
                    'B', zeros(1, 3), ...   % coefficients of the numerator
                    'A', zeros(1, 2));      % coefficients of the denominator
                
switch order
    case 1
        a = 0;
        b = 0;
        c = 4 * bandwidth;
    case 2
        omega0 = bandwidth / 0.53;
        a = 0;
        b = omega0.^2;
        c = 1.414 * omega0;
    case 3
        omega0 = bandwidth / 0.7845;
        a = omega0^3;
        b = 1.1 * omega0.^2;
        c = 2.4 * omega0;
end

LoopFilter.B( 1 ) = a * Tc^2 / 4 + b * Tc /2 + c;
LoopFilter.B( 2 ) = a * Tc^2 / 2             - 2*c;
LoopFilter.B( 3 ) = a * Tc^2 / 4 - b * Tc /2 + c;

LoopFilter.A( 1 ) = -2;
LoopFilter.A( 2 ) = 1;

% Set the initial value
LoopFilter.outputValues( 1 ) = initialValue;
LoopFilter.outputValues( 2 ) = initialValue;