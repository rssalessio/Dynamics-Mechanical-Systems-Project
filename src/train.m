clear all, close all, clc

load('bridgeStructure.txt_mkr.mat');

Mff = M(1:207,1:207); Kff = K(1:207,1:207);

d = 26; %distance between loads
L = 70; %bridge length

V = 34; %speed of the load [m/sec]
%ResonantVelocity = 51.0742; V = ResonantVelocity;


trains = 4; %number of trains/cars

CrossingTime = L/V %time taken to cross the whole bridge
Delay = d/V %time between two cars
MaxTrainOnBridge = L/d %maximum numbers of cars on the bridge
Force = 20e3*9.81;

%% shapes and frequencies
[PHI, natfreq] = eig(Mff\Kff);
natfreq = sqrt(diag(natfreq))/2/pi;

PHI = PHI'; %row <=> i-th shape
SORT = [natfreq, PHI];
SORT = sortrows( SORT, 1 );
natfreq = SORT(:,1);
PHI = SORT(:,2:208);
PHI = PHI';

%principal coordinates
Mq = PHI'*Mff*PHI; Mq = diag( diag(Mq) );
Kq = PHI'*Kff*PHI; Kq = diag( diag(Kq) );

%%

dt = 1e-3; 
TMAX = CrossingTime + (trains-1)*Delay; %smart axis

t=0:dt:TMAX; 
tcross = 0:dt:CrossingTime;

Q = zeros(1, length(t)); %Q is the related to the first mode of vibration

%build loads sequences
for j=1:trains
    Qnext = Force*sin( pi*V*tcross/L );
    
    s(j) = 1 + (j-1)*Delay/dt;      %start time j-th train
    e(j) = s(j)+length(tcross)-1;  %end time j-th train
    
    Q(s(j):e(j)) = Q(s(j):e(j)) + Qnext;
end

figure; plot(t,Q); grid; title('Long Train');

Qkernel = Q(s(3):s(4)-1); %extract period kernel
figure; plot(t(s(3):s(4)-1), Qkernel);

%% fft
N = length(Qkernel);
FFT = fft(Qkernel);

modQ(1) = abs(FFT(1))/N;
modQ(2:N/2) = abs(FFT(2:N/2))*2/N;
phQ = angle(FFT(1:N/2));

DISPLAY = 10;
df = 1/Delay; %armonics frequency
fQ = 0:df:df*(N/2-1);

figure; bar(fQ(1:DISPLAY), modQ(1:DISPLAY), 0.1);

for j=1:DISPLAY
    ResonantVelocity(j) = d*natfreq(1)/j;
end

disp(ResonantVelocity);

%%
IMP_t = zeros(1, length(t));

freq = 0:0.01:10;

for j=1:length(freq)
    omega = 2*pi*freq(j);
    q = ( -omega^2*Mq + Kq )\[1;zeros(206,1)];
    x = PHI*q;
    
    mod(j) = abs( x(idb(8,2)) );
    ph(j) = angle( x(idb(8,2)) );
    
    IMP_t = IMP_t + mod(j)*cos(omega*t + ph(j) );
end

figure; subplot 211; plot(t, IMP_t); title('Impulse Response'); %notice undamped
        subplot 212; plot(freq, mod); title('Spectrum Response'); %only first exc taken into account
   
%%
A_t = zeros(1, length(t));

for j=1:N/2
    omega = 2*pi*fQ(j);
    Q0 = modQ(j)*exp(i*phQ(j));
    q = ( -omega^2*Mq + Kq )\[1;zeros(206,1)]*Q0;
    x = PHI*q;
    
    mod(j) = abs( x(idb(8,2)) );
    ph(j) = angle( x(idb(8,2)) );
    
    A_t = A_t + mod(j)*cos(omega*t + ph(j) );
end

figure; plot(t, A_t); title( strcat('Response to Long Train V=', num2str(V), 'm/s')); %notice undamped
       % subplot 212; bar(fQ(1:DISPLAY), mod(1:DISPLAY), 0.5); title('Spectrum Response'); %only first exc taken into account
