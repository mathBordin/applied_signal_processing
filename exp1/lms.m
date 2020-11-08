%Algoritmo LMS
%Matheus Bordin Gomes
function [W,e,y]=lms(u, d, M, N, mu)
    u_aux=zeros(M,1);
    W=zeros(N,M);
    y=zeros(N,1);
    e=zeros(N,1);
    for n=1:N
        u_aux=[u(n);u_aux(1:M-1)];
        y(n)=W(n,:)*u_aux ;
        e(n)=d(n)-y(n);
        W(n+1,:)=W(n,:)+mu*e(n)*u_aux';
    end
end