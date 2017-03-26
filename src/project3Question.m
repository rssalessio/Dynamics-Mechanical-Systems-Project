function [] = project3Question(nf,NA,NB,Mff,Rff,Kff)

    df = 0.01;
    frequencies = (0:df:15)';
    %----------
    F= zeros(nf,1); F(NA,1) = 1; 
    for i=1:size(frequencies,1);
       omega = frequencies(i)*2*pi;
       x = (-omega^2*Mff+j*omega*Rff+Kff)\F;
       xa = x(NA);
       xb = x(NB);
       moda(i) = abs(xa);
       modb(i) = abs(xb);

       fasa(i) = angle(xa);
       fasb(i) = angle(xb);
    end

    figure;
    subplot 211; plot(frequencies, moda); grid; xlabel('[Hz]'); ylabel('[M]'); title('Y node A - Force on A');
    subplot 212; plot (frequencies, fasa*180/pi); grid; xlabel('[Hz]'); ylabel('[Deg]');

    figure;
    subplot 211; plot(frequencies, modb); grid; xlabel('[Hz]'); ylabel('[M]'); title('Y node B - Force on A');
    subplot 212; plot (frequencies, fasb*180/pi); grid; xlabel('[Hz]'); ylabel('[Deg]');

    %---------
    F= zeros(nf,1); F(NB,1) = 1; 
    for i=1:size(frequencies,1);
       omega = frequencies(i)*2*pi;
       x = (-omega^2*Mff+j*omega*Rff+Kff)\F;
       xa = x(NA)*(-omega^2);
       xb = x(NB)*(-omega^2);
       moda(i) = abs(xa);
       modb(i) = abs(xb);

       fasa(i) = angle(xa);
       fasb(i) = angle(xb);
    end

    figure;
    subplot 211; plot(frequencies, moda); grid; xlabel('[Hz]'); ylabel('[M]'); title('Y acceleration node A - Force on B');
    subplot 212; plot (frequencies, fasa*180/pi); grid; xlabel('[Hz]'); ylabel('[Deg]');

    figure;
    subplot 211; plot(frequencies, modb); grid; xlabel('[Hz]'); ylabel('[M]'); title('Y acceleration node B - Force on B');
    subplot 212; plot (frequencies, fasb*180/pi); grid; xlabel('[Hz]'); ylabel('[Deg]');
end