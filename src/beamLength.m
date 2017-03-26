function [L] = beamLength(parameters, frequency, approximationType, approximationParam)
%
% Parameters is a vector with parameters P,A,J,E -> parameters = [P,A,J,E];
% Frequency denotes the frequency range to be approximated
% Approximation is a scale factor that multiplied by the frequency gives
% the approximation, by default is 1 decade (10)
%
% The formula is 2*pi*f = (pi/L)^2   * sqrt(E*J/M);
% Therefore L = sqrt(pi/(2*f)  * sqrt(E*J/M));

% switch(nargin)
%     case 2
%         approximation =10 ;
%     case 3 ;
%     otherwise
%         error('beamLength:TooManyInputs', ...
%         'requires at least 2 inputs');
% end

frequency = abs(frequency);
parameters = abs(parameters);
switch(approximationType)
    case ApproxType.HalfPower
        fmax = frequency*sqrt((1+sqrt(2)));
    case ApproxType.FreqRange
        fmax = frequency*approximationParam;
    case ApproxType.DerivativeRule
        a=frequency;
        k=approximationParam;
        r = roots([-k,0,3*a^2*k,0,-3*a^2*k,-4*a^2,k*a^6]);
        r = round(r*10000)/10000.0;
        fmax = r(real(r)>a & imag(r)==0);        
    otherwise
        fmax = frequency*10;
end


 
M = parameters(1)*parameters(2);
sq1 = sqrt(parameters(4)*parameters(3)/M);
L = sqrt(pi*sq1/(2*fmax));


end