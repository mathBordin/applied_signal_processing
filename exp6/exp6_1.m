%%
% Experiência 6: Beamforming
% Matheus Bordin Gomes
clear;clc;close all;

%%
%Dados iniciais
c = 3*10^8;         %Velocidade da luz no vácuo
f = 500*10^9;        %Frequência do sinal recebido
lambda = c/f;       %Comprimento de onda do sinal recebido
d = lambda/2;       %Distância entre as antenas
theta_0 = 30;        %Ângulo de incidência do sinal
theta = -90:1:90;   %Variação do ângulo de incidência
M = 5;              %Número de antenas

%%
%Amostragem
fa = 2*f;           %Frequência de amostragem
Ta = 1/fa;          %Período de amostragem

%%
%B(teta,teta_0)
omega = 2*pi*f*Ta;
u_theta = d*sind(theta)/(Ta*c);
u_theta_0 = d*sind(theta_0)/(Ta*c);
B=(1/M)*(1-exp(1i*omega*M*(u_theta-u_theta_0)))./(1-exp(1i*omega*(u_theta-u_theta_0)));

%Gráfico de |B(teta,teta_0)|
figure(1);
subplot(2,2,1);
plot(theta,abs(B));
title('|B(\theta,20º)| para d=\lambda/4');
xlabel('\theta');

%%
%Gráficos de B(teta,teta_0) para outros valores de d

%%
%d = lambda
d = lambda;
u_theta = d*sind(theta)/(Ta*c);
u_theta_0 = d*sind(theta_0)/(Ta*c);
B=(1/M)*(1-exp(1i*omega*M*(u_theta-u_theta_0)))./(1-exp(1i*omega*(u_theta-u_theta_0)));

subplot(2,2,2);
plot(theta,abs(B));
title('|B(\theta,20º)| para d=\lambda');
xlabel('\theta');

%%
%d = 3lambda/4
d = 3*lambda/4;
u_theta = d*sind(theta)/(Ta*c);
u_theta_0 = d*sind(theta_0)/(Ta*c);
B=(1/M)*(1-exp(1i*omega*M*(u_theta-u_theta_0)))./(1-exp(1i*omega*(u_theta-u_theta_0)));

subplot(2,2,3);
plot(theta,abs(B));
title('|B(\theta,20º)| para d=3\lambda/4');
xlabel('\theta');

%%
%d = lambda/2
d = lambda/2;
u_theta = d*sind(theta)/(Ta*c);
u_theta_0 = d*sind(theta_0)/(Ta*c);
B=(1/M)*(1-exp(1i*omega*M*(u_theta-u_theta_0)))./(1-exp(1i*omega*(u_theta-u_theta_0)));

subplot(2,2,4);
plot(theta,abs(B));
title('|B(\theta,20º)| para d=\lambda/2');
xlabel('\theta');

%%
%Interferência de outro sinal
A1 = 1;
A2 = 0.5;
theta_2 = theta;
amp_inter = abs(A1+B*A2);
figure(2);
plot(theta_2,amp_inter);
title('Amplitude do sinal incidente para d=\lambda/2');
xlabel('\theta_2');