%Experiência 1 - Filtros Adaptativos
%Matheus Bordin Gomes

clear; close all; clc; 

%%%
%Equalização Adaptativa
%%%

q_v = [2.9 3.1 3.3 3.5];
MSE_lms_estim = zeros(4,1);
MSE_lms_estim_dB = zeros(4,1);
MSE_nlms_estim = zeros(4,1);
MSE_nlms_estim_dB = zeros(4,1);
for p = 1:4
    N = 500; %Número de amostras do sinal
    n_amostras = 1:1:N;
    M = 11; %Número de coeficientes do equalizador adaptativo
    N_exp = 100;
    H_M = 3;

    a = 2*round(rand(1,N))-1; %Sina binário i.i.d. transmitido
    nd = 6;

    k = 1:H_M;
    q = q_v(p);
    h = 0.5*(1+cos(2*pi*(k-2)/q)); %Resposta do canal

    var_n = 0.001;
    n = sqrt(var_n)*randn(1,N); %AWGN no sinal de saída do canal

    u = filter(h,1,a) + n; %Sinal recebido
    delay_zeros = zeros(1,nd);
    d = [delay_zeros a(1:end-nd)];

    %Cálculo da matriz de autocorrelação do sinal u
    ru_0 = h(1)^2 + h(2)^2 + h(3)^2;
    ru_1 = h(1)*h(2) + h(2)*h(3);
    ru_2 = h(1)*h(3);

    R = zeros(M);
    for i = 1:M
        for j = 1:M
            if i==j
                R(i,j) = ru_0;
            elseif (abs(i-j) == 1)
                R(i,j) = ru_1;
            elseif (abs(i-j) == 2)
                R(i,j) = ru_2;                
            end
        end
    end

    %Autovalores da matriz R de interesse
    eigenvalues = eig(R); 
    eig_max = max(eigenvalues);
    eig_min = min(eigenvalues);
    eig_max_min = eig_max/eig_min;

    %Algoritmo LMS
    mu_lms = 0.075;
    MSE_lms = zeros(N,1);
    for i = 1:100
        [W,e,y]=lms(u, d, M, N, mu_lms);
        MSE_lms = MSE_lms + e.^2;    
    end
    MSE_lms = MSE_lms/N_exp;
    MSE_lms_estim(p) = mean(MSE_lms(0.95*N:end));
    MSE_lms_estim_dB(p) = 10*log10(mean(MSE_lms(0.95*N:end)));
    
    figure(1);
    if p == 1
        plot(n_amostras,10*log10(MSE_lms),'b');
    elseif p == 2
        plot(n_amostras,10*log10(MSE_lms),'r');
    elseif p == 3
        plot(n_amostras,10*log10(MSE_lms),'g');
    elseif p == 4
        plot(n_amostras,10*log10(MSE_lms),'y');
    end
    xlabel('Amostras');
    ylabel('EMS(dB)');
    title('LMS');
    hold on;

    %Algoritmo NLMS
    mu_nlms = 0.7; 

    del = 10^(-5);
    MSE_nlms = zeros(N,1);
    for i = 1:100
        [W,e,y]=nlms(u, d, M, N, mu_nlms,del);
        MSE_nlms = MSE_nlms + e.^2;    
    end
    MSE_nlms = MSE_nlms/N_exp;
    MSE_nlms_estim(p) = mean(MSE_nlms(0.95*N:end))
    MSE_nlms_estim_dB(p) = 10*log10(mean(MSE_nlms(0.95*N:end)));
    
    figure(2);
    if p == 1
        plot(n_amostras,10*log10(MSE_nlms),'b');
    elseif p == 2
        plot(n_amostras,10*log10(MSE_nlms),'r');
    elseif p == 3
        plot(n_amostras,10*log10(MSE_nlms),'g');
    elseif p == 4
        plot(n_amostras,10*log10(MSE_nlms),'y');
    end
    xlabel('Amostras');
    ylabel('EMS(dB)');
    title('NLMS');
    hold on;
    
    if p == 1
        figure(3);        
        plot(n_amostras,10*log10(MSE_lms),'b');
        xlabel('Amostras');
        ylabel('EMS(dB)');
        title('LMS');
        
        figure(4);
        plot(n_amostras,10*log10(MSE_nlms),'b');
        xlabel('Amostras');
        ylabel('EMS(dB)');
        title('NLMS');
    end
end

MSE_lms_estim
MSE_lms_estim_dB
MSE_nlms_estim 
MSE_nlms_estim_dB 

figure(1);
legend('q = 2.9', 'q = 3.1', 'q = 3.3', 'q = 3.5');
figure(2);
legend('q = 2.9', 'q = 3.1', 'q = 3.3', 'q = 3.5');
