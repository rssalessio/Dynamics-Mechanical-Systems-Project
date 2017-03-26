function []= project5Question2(sysomega,nf,nc,M,R,K,alpha,beta)
    [phi,~]=eig(M\K);
    R = alpha*M+beta*K;
    M = phi'*M*phi;
    K = phi'*K*phi;
    R = phi'*R*phi;
    
end