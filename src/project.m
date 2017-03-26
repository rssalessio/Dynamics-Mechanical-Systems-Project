clc; clear all;
load('bridgeStructure.txt_mkr.mat');
%load('bridgeStructure.New_mkr.mat');
%load('bridgeStructure.Newb_mkr.mat');
%load('bridgeStructure.springs_mkr.mat');
%load('bridgeNew.txt_mkr.mat');
%load('bridgetestlast.txt_mkr.mat');
close all;
% Node O1: 1
% Node A: 32
% Node B: 13
% Node O2: 64
alpha=0.2;
beta=1e-4;
A = 32;%50; %32
B = 13;%20;%13;
O1 = 1;
O2 =64;%101;%64;

NA = idb(A,2);
NB = idb(B,2);
N01 = idb(O1,2);
N02 = idb(O2,2);

n = size(idb,1)*3;
nc = 3;
nf = n-nc;

Mff = M(1:nf,1:nf);
Kff = K(1:nf,1:nf);
Rff = R(1:nf,1:nf);

Mfc = M(1:nf, nf+1:nf+nc);
Kfc = K(1:nf, nf+1:nf+nc);
Rfc = R(1:nf, nf+1:nf+nc);

Mcf = M(nf+1:nf+nc, 1:nf);
Kcf = K(nf+1:nf+nc, 1:nf);
Rcf = R(nf+1:nf+nc, 1:nf);
Mcc = M(nf+1:nf+nc, nf+1:nf+nc);
Kcc = K(nf+1:nf+nc, nf+1:nf+nc);
Rcc = R(nf+1:nf+nc, nf+1:nf+nc);
%%
[v,d]=eig(Mff\Kff);
[freq]=sqrt(diag(d))/2/pi;
v = v'; %row <=> i-th shape
SORT = [freq, v];
SORT = sortrows( SORT, 1 ); 
freq = SORT(:,1);
v = SORT(:,2:nf+1);
v = v';


% Anew = [zeros(nf,nf), eye(nf,nf); -Mff\Kff, -Mff\Rff];
%  [shapes2,puls2]=eig(Anew);
%  
%  [omega,I]=sort( diag(puls2) );
%  freq = unique(abs(imag(omega)/(2*pi)));
%%  
%project3Question(nf,NA,NB,Mff,Rff,Kff);
%%
%project4Question(Mff,Mfc,Rff,Rfc,Kff,Kfc,NA,NB);
%%
% 
% [X,Y]=meshgrid(1:100, 1:1:20);
% freq=freq';
%  surf(26*freq(Y)./(2*pi*X));
%  xlabel('k'); ylabel('freq'); zlabel('V');
% k=1:1:10;
% i=6;
% %v = 26*freq(i)./(k);
% %plot(k,v);
% v=1:0.0001:90;
% v=1./v;
% freqk=freq(2);
% Cn = 26*v'*freqk';
% kmax = ceil(max(Cn(:)));
%  sCn=size(Cn);
%  kvals = [0];
%  vvals = [0];
% for k=1:kmax
%    eps = abs(ones(sCn)*k - Cn);
%    if (min(eps(:)) < 1e-5 *13*45/pi)
%        [row,col] = find(eps == min(eps(:)));
%         kvals = [kvals , k];
%         vvals = [vvals , 1/v(row)];
%        disp(['K: ' num2str(k) '- Matches: ' num2str(length(find(eps<1e-5*13*45/pi))) ' - [MIN]  Diff: ' num2str(eps(row,col)) ' - Velocity: ' num2str(1/v(row)) ' - Frequency: ' num2str(freq(col))]);
%    end
% end
% 
% i=2;
% Vmax = 89.4965;
% k=1:1:max(kvals);
% y = Vmax./(1+(k-i)*i^(-1));
% k=[k];
% y=[y];
% figure;
% plot(kvals,vvals,'o--r', k,y); legend('Resonance velocities','Resonance velocities from law'); title('2nd frequency velocities'); xlabel('k'); ylabel('V');
% axis([0,60, 0,Vmax]);
% 
% % y2 = 51./xx;

V=14;
d=26;
L=70;
ntrain=4;

T=L/V;
tau=d/V;
TMAX = T+ (ntrain-1)*tau;
t=0:0.1:TMAX;
figure;
y=zeros(size(t));
for i=1:ntrain
   Q = sin(pi*V*(t-(i-1)*tau)/L);
   Q(t>(T+i*tau)) = 0;
   Q(t< (i-1)*tau) = 0;
   Q(find(Q<0))=0;
   plot(t,Q,'--'); hold on; xlabel('t'); ylabel('y');  hold on;
   y=y+Q;
end
plot(t,y);

N = length(y);
    
    %--- Task 5: FFT ---%
df = 1/tau;
fmax= df*(N/2-1);
input_nome_frequencies = 0:df:fmax;
fftout = fft(y);
fftabs = zeros(1,N/2);
fftang = zeros(1,N/2);

fftabs(1) = abs(fftout(1))/N;
fftabs(2:N/2) = abs(fftout(2:N/2))*2/N;
fftang(1:N/2) = angle(fftout(1:N/2));

figure;
    subplot 211; bar(input_nome_frequencies, fftabs,0.1); xlabel('[Hz]'); ylabel('[ ]'); title('FFT Input Signal: ');
    subplot 212; bar(input_nome_frequencies, fftang*180/pi,0.1); xlabel('[Hz]'); ylabel('[Deg]');
    
 Mff = v'*Mff*v; Mff=diag(Mff);
 Kff = v'*Kff*v; Kff = diag(Kff);
 
IMP_t = zeros(1, length(t));

freq = 0:0.01:10;
A_t = zeros(1, length(t));

for j=1:N/2
    omega = 2*pi*input_nome_frequencies(j);
    Q0 = fftabs(j)*exp(i*fftang(j));
    q = ( -omega^2*Mff+ Kff )\[1;zeros(nf-1,1)]*Q0;
    x = v*q;
    
    mod(j) = abs( x(30) );
    ph(j) = angle( x(30) );
    
    A_t = A_t + mod(j)*cos(omega*t + ph(j) );
end

figure; subplot 211; plot(t, A_t); title( strcat('Response to Long Train V=', num2str(V), 'm/s')); %notice undamped
        %subplot 212; bar(fQ(1:DISPLAY), mod(1:DISPLAY), 0.5); title('Spectrum Response'); %only first exc taken into account

% figure;
% plot(kvals,vvals,'o--b');hold on;
% scatter(xx,y2,'r');
% vvals(2)
% 
% mv=[51.223,89.4965,106.0798,186.52,186.6,193.7,198.37,208.13];
% figure;plot(mv);hold on; plot(freq(1:8));
% 
% figure;
% mv =[mv(2)/mv(1), mv(3)/mv(2), mv(4)/mv(3), mv(5)/mv(4),mv(6)/mv(5),mv(7)/mv(6),mv(8)/mv(7)];
% mf = [freq(2)/freq(1),freq(3)/freq(2),freq(4)/freq(3),freq(5)/freq(4),freq(6)/freq(5),freq(7)/freq(6),freq(8)/freq(7)];
% plot(mv); hold on; plot(mf);
% 
% 
% [x,resnorm] = lsqcurvefit(@myfun,[1;1;1;1;1],mv,mf);
% 
% k=0:0.1:2;
% y = x(1)./mv +x(2)./mv.^2;
% figure;
% plot(y); hold on; plot(mf);
% 

%%
% close all;
% clear all;
% V=3;
% L=70;
% d=26;
% T=L/V;
% tau=d/V;
% ntrain=2;
% figure;
% t = 0:0.1:100;
% Q = zeros(size(t));
% for i=0:ntrain-1
%     y = sin(pi*V*t/L-i*8.6);%y(t>T+i*tau)=0; y(t<i*tau)=0;
%     plot(t,y); hold on;
% end

%  for k=1:10
%     26*freq(1)/(2*pi*k)  
%  end

%project5Question(sysomega,nf,nc,Mff,Mfc,Mcf,Mcc,Rff,Rfc,Rcf,Rcc,Kff,Kfc,Kcf,Kcc);
%project5Question2(sysomega,nf,nc,M,R,K,alpha,beta);
%project6Question(sysomega,nf,nc,Mff,Mfc,Rff,Rfc,Kff,Kfc,alpha,beta
