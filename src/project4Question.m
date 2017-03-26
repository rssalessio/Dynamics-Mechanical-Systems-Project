function [] = project4Question(Mff,Mfc,Rff,Rfc,Kff,Kfc,NA,NB)
    quake=load('sismaspost.txt');
    plotQuake(quake);
    t= (quake(:,1))';
    yO1 = quake(:,2);
    yO2 = quake(:,3);

    N = size(quake,1);
    T = quake(end,1);
    fftO1 = fft(yO1);
    fftO2 = fft(yO2);

    df=1/T;
    fmax = df*(N/2-1);

    frequencies = (0:df:fmax)';

    fftabsO1(1) = abs(fftO1(1))/N;
    fftabsO2(1) = abs(fftO2(1))/N;
    fftabsO1(2:N/2) = abs(fftO1(2:N/2))*2/N;
    fftabsO2(2:N/2) = abs(fftO2(2:N/2))*2/N;
    fftfasO1(1:N/2) = angle(fftO1(1:N/2));
    fftfasO2(1:N/2) = angle(fftO2(1:N/2));

    figure;
    subplot 211; bar(frequencies, fftabsO1); xlabel('[Hz]'); title('FFT Time Series O1 - Quake');  axis([-2 10 0 0.03]);
    subplot 212; bar(frequencies, fftfasO1*180/pi); xlabel('[Hz]'); ylabel('[Deg]');  axis([-2 10 -200 200]);


    figure;
    subplot 211; bar(frequencies, fftabsO2); xlabel('[Hz]'); title('FFT Time Series O2 - Quake'); axis([-2 10 0 0.03]);
    subplot 212; bar(frequencies, fftfasO2*180/pi); xlabel('[Hz]'); ylabel('[Deg]');  axis([-2 10 -200 200]);

    qrO1 = zeros(1,N);
    yA = zeros(1,N);
    yB = zeros(1,N);
    yApp = zeros(1,N);
    yBpp = zeros(1,N);
    for i =1:length(frequencies)
        omega = frequencies(i)*2*pi;
        XcO1 = fftabsO1(i)*exp(j*fftfasO1(i));
        XcO2 = fftabsO2(i)*exp(j*fftfasO2(i));
        Xc = [0; XcO1; XcO2;0];
        Qc = (-omega^2*Mfc+j*omega*Rfc+Kfc)*Xc;
        x = (-Mff*omega^2 +j*omega*Rff+Kff)\(-Qc);

        xa = x(NA);
        xb = x(NB);
        xapp = x(NA)*(-omega^2);
        xbpp = x(NB)*(-omega^2);

        moda(i) = abs(xa);
        modb(i) = abs(xb);
        modapp(i) = abs(xapp);
        modbpp(i) = abs(xbpp);


        fasa(i) = angle(xa);
        fasb(i) = angle(xb);
        fasapp(i) = angle(xapp);
        fasbpp(i) = angle(xbpp);

        qrO1 = qrO1+fftabsO1(i)*cos(omega*t+fftfasO1(i));
      %  if (frequencies(i) < 0.5)
            yA = yA+ moda(i)*cos(omega*t+fasa(i));
            yApp = yApp+ modapp(i)*cos(omega*t+fasapp(i));
            yB = yB+ modb(i)*cos(omega*t+fasb(i));
            yBpp = yBpp+ modbpp(i)*cos(omega*t+fasbpp(i));
       % end
    end
    
    figure;
    plot(t, qrO1); grid; xlabel('Time [s]'); ylabel('[M]'); title('Quake in O1 rebuilt');

    figure;
    subplot 211; plot(frequencies, moda); grid; xlabel('[Hz]'); title('Spectrum yA - Quake'); axis([-2 10 0 0.03]);
    subplot 212; plot(frequencies, fasa*180/pi); grid; xlabel('[Hz]'); ylabel('[Deg]');  axis([-2 10 -200 200]);


    figure;
    subplot 211; plot(frequencies, modapp); grid; xlabel('[Hz]'); title('Spectrum yApp - Quake'); axis([-2 10 0 2]);
    subplot 212; plot(frequencies, fasapp*180/pi); grid; xlabel('[Hz]'); ylabel('[Deg]');   axis([-2 10 -200 200]);

    figure;


    subplot 211; plot(frequencies, modb); grid; xlabel('[Hz]'); title('Spectrum yB - Quake'); axis([-2 10 0 0.03]);
    subplot 212; plot(frequencies, fasb*180/pi); grid; xlabel('[Hz]'); ylabel('[Deg]');  axis([-2 10 -200 200]);

    figure;


    subplot 211; plot(frequencies, modbpp); grid; xlabel('[Hz]'); title('Spectrum yBpp - Quake'); axis([-2 10 0 2]);
    subplot 212; plot(frequencies, fasbpp*180/pi); grid; xlabel('[Hz]'); ylabel('[Deg]');  axis([-2 10 -200 200]);

    figure;
    subplot 211; plot(t, yA);grid;xlabel('Time [s]'); ylabel('[M]'); title('yA - Quake');
    subplot 212; plot(t,yApp);grid;xlabel('Time [s]'); ylabel('[M]'); title('yApp - Quake');

    figure;
    subplot 211; plot(t, yB);grid;xlabel('Time [s]'); ylabel('[M]'); title('yB - Quake');
    subplot 212; plot(t,yBpp);grid;xlabel('Time [s]'); ylabel('[M]'); title('yBpp - Quake');

    figure;
    subplot 211; plot(t,abs(yA-yB)); grid;xlabel('Time [s]'); ylabel('[M]'); title('abs(yA-yB) - Quake');
    subplot 212; plot(t,abs(yApp-yBpp)); grid;xlabel('Time [s]'); ylabel('[M]'); title('abs(yApp-yBpp) - Quake');
    disp(['Maximum yA: ' num2str(max(yA))]);
end