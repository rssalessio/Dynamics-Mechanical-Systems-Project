function []= project6Question(sysomega,nf,nc,Mff,Mfc,Rff,Rfc,Kff,Kfc,alpha,beta)
    [phi,~]=eig(Mff\Kff);
    Rff = alpha*Mff+beta*Kff;
    Mff = phi'*Mff*phi;
    Kff = phi'*Kff*phi;
    Rff = phi'*Rff*phi;
    i=146;
    %2hz-> i=146
    disp(['Mff componente 2 hertz: ' num2str(Mff(i,i)) ' - Kff: ' num2str(Kff(i,i))]);
    
end