function [] = project5Question(sysomega,nf,nc,Mff,Mfc,Mcf,Mcc,Rff,Rfc,Rcf,Rcc,Kff,Kfc,Kcf,Kcc)
    d = 26;
    M = 1000;
    g = 9.81;
    V = 1:1:150;
    V = V'/3.6;

    InvV= 1./V;
    C = 13*InvV*sysomega'/pi;

    Nv = length(V);
    Q = zeros(nf,1);
    Fy = zeros(nc,1);
    t=0:0.1:100;
    t=t';
    for i = 150:150
        yO1 = zeros(length(t),1);
        FO1 = zeros(length(t),1);
        for k=1:1000  
           omega = k*2*pi*V(i)/d;
           [res,I] = max(abs (omega*ones(size(sysomega))-sysomega));
           if (res <20)
              disp([ num2str(omega) '- ' num2str(sysomega(I))]); 
           end


           A1 = (-omega^2*Mff+j*omega*Rff+Kff);
           B1 = (-omega^2*Mfc+j*omega*Rfc+Kfc);
           A2 = (-omega^2*Mcf +j*omega*Rcf+Kcf);
           B2 = (-omega^2*Mcc+j*omega*Rcc+Kcc);
           Fy = [0;1;0;0]*M*g*V(i)/d;
           A = [A1, B1; A2, B2];
           x = A\[Q;Fy];
           xc = x(nf+1:end);

           FO1 = FO1+M*g*cos(omega*t)*V(i)/d;
           yO1 = yO1+abs(xc(2))*cos(omega*t+angle(xc(2)));
        end
        figure;
       subplot 211; plot(t,FO1); axis; xlabel('Time [s]'); ylabel('[N]');
        subplot 212; plot(t,yO1); axis; xlabel('Time [s]'); ylabel('[M]'); hold on;
    end 
end