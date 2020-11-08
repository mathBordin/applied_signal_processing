%Algoritmo NLMS
%Matheus Bordin Gomes
function [W,e,y]=nlms(u, d, M, N, mu_0, delta) 
    u_aux=zeros(M,1);
    mu = 0;
    
    W=zeros(N,M);
    y=zeros(N,1);
    e=zeros(N,1);
    
    for n=1:N
        u_aux=[u(n);u_aux(1:M-1)];
        y(n)=W(n,:)*u_aux;
        e(n)=d(n)-y(n);
        W(n+1,:)=W(n,:)+mu*e(n)*u_aux';
        mu = mu_0/(delta + norm(u_aux)^2);
    end
end

