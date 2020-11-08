%%
% Experiência 6: Beamforming
% Matheus Bordin Gomes
clear;clc;close all;

%%
%Dados iniciais
c = 3*10^8;         %Velocidade da luz no vácuo
f = 60*10^9;        %Frequência do sinal recebido
Omega = 2*pi*f;     %Frequência angular do sinal recebido
lambda = c/f;       %Comprimento de onda do sinal recebido
d = lambda/2;       %Distância entre as antenas
M = 8;              %Número de antenas
theta_0 = 20;       %Ângulo de incidência do sinal de interesse
load('sinais.mat'); %Sinal recebido em cada uma das oito antenas

%%
% Ruído adicional para o item 6. Caso não queira ruído, deixar noise_amp em zero. 
noise_amp = 10;  
noise = noise_amp*randn(size(sinal_recebido));
sinal_recebido = sinal_recebido + noise;

%%
%Amostragem
fa = 10^12;           %Frequência de amostragem
Ta = 1/fa;          %Período de amostragem

%%
%Projeto do filtro Butterworth passa-baixa
[b,a] = butter(6,f/(fa/2));
figure(1);
freqz(b,a,1024,fa);
title('Resposta em frequência do filtro Butterworth passa-baixas');

%%
%Decodificação dos sinais xim e xqm em cada antena
arg_t = ((1:size(sinal_recebido,1))-1)/fa;
%xi recebido na antena 1
xim = filter(b,a,2*cos(Omega*arg_t)'.*sinal_recebido);
figure(2);
subplot(2,1,1);
plot(xim(:,1));
title('xi_1');
xlabel('Amostras');
ylabel('Amplitude');

%xq recebido na antena 1
xqm = filter(b,a,-2*sin(Omega*arg_t)'.*sinal_recebido);
subplot(2,1,2);
plot(xqm(:,1));
title('xq_1');
xlabel('Amostras');
ylabel('Amplitude');

xm = xim + 1j*xqm;

%%
%Projeto dos coeficientes do vetor w
m = 0:1:(M-1);
v_h = conj(exp(1j*m*Omega*d*sind(theta_0)/c));
w_h = (1/M)*v_h;

%%
%Obtenção do sinal recebido
y = w_h*xm.';
%Parte real de y
figure(3);
subplot(2,2,1);
plot(real(y));
title('R(y)');
xlabel('Amostras');
ylabel('Amplitude');
%Parte imaginária de y
subplot(2,2,2);
plot(imag(y));
title('Im(y)');
xlabel('Amostras');
ylabel('Amplitude');

%%
%Recuperação dos símbolos transmitidos
N = round(10^(-9)*fa);
n_simbols = round(length(y(1,:))/N);
y_rcv=zeros(n_simbols,1);
for k = 1:n_simbols
    y_rcv(k) = mean(y(1+(k-1)*N:k*N));
end

%Parte real de y_rcv
y_rcv_plot = y_rcv;
y_rcv_plot(k+1)=0;
x_plot = 0:10;
subplot(2,2,3);
stairs(x_plot,real(y_rcv_plot));
title('R(y) - media');
xlabel('Amostras');
ylabel('Amplitude');
%Parte imaginária de y_rcv
subplot(2,2,4);
stairs(x_plot,imag(y_rcv_plot));
title('Im(y) - media');
xlabel('Amostras');
ylabel('Amplitude');

%%
%Decodificação dos símbolos
symbols_16_qam = qammod(0:15,16);
[~,ind] = min(abs(symbols_16_qam - y_rcv),[],2);
simbols = symbols_16_qam(ind);

%%
% Erros
load('mensagem.mat');
message = symbols_16_qam(mensagem+1);

erros = sum(sum(simbols ~= message))  %#ok<NOPTS>

figure(4);
scatter(real(message), imag(message));
hold on;
scatter(real(simbols), imag(simbols), 'x');
hold off;
title('Scatter dos símbolos tx e rx');
legend('Mensagem transmitida', 'Mensagem recebida');