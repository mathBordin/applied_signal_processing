%%
% Prova Final - PSI3531 (Processamento de Sinais Aplicado)
% Questão 1. Versão com ângulos com erros.
% Matheus Bordin Gomes

clear;clc;close all;

%%
%Dados iniciais
c = 3*10^8;             %Velocidade da luz no vácuo
f = 60*10^9;            %Frequência do sinal recebido
Omega = 2*pi*f;         %Frequência angular do sinal recebido
lambda = c/f;           %Comprimento de onda do sinal recebido
d = lambda/2;           %Distância entre as antenas
M = 8;                  %Número de antenas
theta = [20 45 -15];    %Ângulo de incidência dos sinais
load('sinais.mat');     %Sinal recebido em cada uma das oito antenas
var_n = 900;            %Variância do ruído nas antenas
fi = [3 ; 1 ; 1];       %Amplitude dos sinais recebidos

%%
% Ruído adicional para o item 6. Caso não queira ruído, deixar noise_amp em zero. 
noise_amp = 0;  
noise = noise_amp*randn(size(sinal_recebido));
sinal_recebido = sinal_recebido + noise;

%%
%Amostragem
fa = 10^12;         %Frequência de amostragem
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
m = (0:1:(M-1))';

%Monta a matriz C, com os vetores v(omega,theta+err)
C = (exp(1j*m*Omega*d*sind(theta + normrnd(0,0.01, size(theta)))/c));

%Monta o vetor V, com os vetores v(omega,theta)
phi = theta ; %Adiciona variação no ângulo de incidência dos sinais
V = (exp(1j*m*Omega*d*sind(phi)/c));

%Calcula sigma
sig = fi*fi';

%Calcula Sr
Sr = var_n*eye(size(V,1));

%Monta g
g_h = zeros(3,1);
g_h(1) = 1;
g = g_h;

%Monta C
w = inv(V*(sig')*V' + Sr)*C*inv(C'*inv(V*(sig')*V' + Sr)*C)*g; %#ok<*MINV>
w_h = w';

%Plota a resposta B
theta_x = -90:0.1:90;
V_theta = (exp(1j*m*Omega*d*sind(theta_x)/c));
B = w_h*V_theta;
figure(3);
plot(theta_x,abs(B));
title('|B(\theta)|');
xlabel('\theta');
ylabel('|B|');
xticks(-90:5:90);
grid;
%%
%Obtenção do sinal recebido
y = w_h*xm.';
%Parte real de y
figure(4);
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

figure(5);
scatter(real(message), imag(message));
hold on;
scatter(real(simbols), imag(simbols), 'x');
hold off;
title('Scatter dos símbolos tx e rx');
legend('Mensagem transmitida', 'Mensagem recebida');