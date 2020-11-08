%Experi�ncia 1 - Filtros Adaptativos
%Matheus Bordin Gomes

clear; close all; 

%%%
%Elimina��o de Interfer�ncias
%%%

N = 500;
n = 1:N;
M = 2;
n_exp = 500;
J_exp1 = zeros(N,3);

%r=xcorr(u,1,'biased') %C�lcula a autocorrela��o experimentalmente
R = [12.5 10.1127 ; 10.1127 12.5]; %Matriz de autocorrela��o de u
p = [2.1651 ; 1.0168]; %Vetor de correla��o cruzada entre d e u
wo = R\p; %Coeficientes �timos do filtro
Jmin = 0.01; %Erro quadr�tico m�nimo

mu_max = 2/max(eig(R)); %Passo m�ximo para o Steepest Descent
mu_max_exp = 0.069; %Passo m�ximo experimental

mu_v = [0.01 0.03 0.05];

%for t = 1:3
    %Roda o algoritimo n_exp vezes 
    t = 2;
    mu = mu_v(t);
    for k = 1:n_exp
        var_s = 0.01;
        s = sqrt(var_s)*randn(1,N); %Sinal que se deseja medir

        phi_v = 2*pi*rand;
        x = sin(2*pi*n/10 + pi/6 + phi_v); %Interfer�ncia a ser eliminar

        phi_u = phi_v;
        u = 5*sin(2*pi*n/10+phi_u); %Sinal correlacionado com a interfer�ncia

        d = s + x; %Resposta desejada do filtro

        w = zeros(M,1);
        w1 = zeros(N,1);
        w2 = zeros(N,1);
        y = zeros(N,1);
        e = zeros(N,1);
        eq = zeros(N,1);

        u_v = [u(1) ; 0];
        for i = 1:N
            y(i) = dot(transpose(u_v),w); %Estima��o do filtro
            e(i) = d(i) - y(i); %Erro de estima��o
            eq(i) =  e(i)^2;
            J_exp1(i,t) = J_exp1(i,t) + eq(i);
            w = w + mu*e(i)*u_v; %Ajusta coeficientes
            w1(i) = w(1);
            w2(i) = w(2);
            if i~=N
                u_v = [u(i+1); u_v(1:M-1, 1)]; %Atualiza valor de u_v
            end
        end
    end

    %Plota sinal correlacionado, erro de estima��o e sinal de interesse
    figure(1);
    plot(n, u, 'r', n, e, 'b', n, s, 'g');
    xlabel('Numero de itera��es');
    legend('Sinal correlacionado','Erro de estima��o', 'Sinal de interesse');

    %Plota coeficientes do filtro
    figure(2);
    plot(n, wo(1)*ones(size(n)), 'r', n, wo(2)*ones(size(n)), 'b', n, w1, 'g', n, w2, 'y');
    xlabel('Numero de itera��es');
    ylabel('Valor dos coeficientes do filtro adaptativo');
    legend('w1 �timo','w2 �timo', 'w1', 'w2');

    %Tra�a as curvas de n�vel da superf�ce de erro e a trajet�ria dos
    %coeficientes
    figure(3);
    clf;
    J = zeros(N);
    w_t1 = linspace(-0.6,0.8,N);
    w_t2 = linspace(-0.6,0.8,N);
    for i = 1:N
        for j = 1:N
            J(i,j)= var(d) - 2*[w_t1(i) w_t2(j)]*p + + [w_t1(i) w_t2(j)]*R*[w_t1(i) ; w_t2(j)];
        end
    end
    contour(w_t1,w_t2,J, linspace(0,0.4,20));
    hold on;
    plot(w2,w1,'b');

    %Plota erro quadr�tico em dB
    figure(4);
    plot(n, 10*log10(eq), 'b');
    xlabel('Numero de itera��es');
    ylabel('Erro quadr�tico (dB)');

    %Estimativa do erro
    J_exp1(:,t) = J_exp1(:,t)/n_exp;
    MSE = mean(J_exp1(0.9*N:end));
    EMSE = MSE - Jmin;
    desajuste = EMSE/Jmin;

    %Plota J(n)
    figure(5);
    plot(n, 10*log10(J_exp1(:,t)), 'b');
    xlabel('Numero de itera��es');
    ylabel('Estimativa de J(n)');
%end

%Plota todos os J(n)
%figure(6);
%plot(n, 10*log10(J_exp1(:,1)), 'b', n, 10*log10(J_exp1(:,2)), 'r', n, 10*log10(J_exp1(:,3)), 'g');
%xlabel('Numero de itera��es');
%ylabel('Estimativa de J(n)');
%legend('mu = 0.01', 'mu = 0.03', 'mu = 0.05');