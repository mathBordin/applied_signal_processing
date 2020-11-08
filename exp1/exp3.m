%Experiência 1 - Filtros Adaptativos
%Matheus Bordin Gomes

clear; close all; clc; 

%%%
%Cancelamento de Eco Ac ?ustico
%%%

[locutor,Fs]=audioread('locutor.wav'); %Carrega sinal de áudio
load('ri.mat'); %Carrega resposta ao impulso do ambiente que gera eco
n=1:1:size(locutor); %Vetor com o número das amostras

%Gera eco do sinal de áudio
eco = filter(ri256, 1, locutor); 
var_noise = 10^(-4);
%interf = sqrt(var_noise)*randn(size(eco)); 
[interf, Fs_2] = audioread('eng.wav'); %Para Double Talk
eco = eco + interf;

%Plota sinais x(n) e d(n)
figure(1)
plot(n,locutor,'b',n,eco,'r');
xlabel('Amostras');
ylabel('Amplitude');
legend('x(n)','d(n)');

%NLMS
mi = 0.1;
delta = 10^(-5);
M = 256;
[W,e,y]=nlms(locutor, eco, M, size(locutor,1), mi, delta);

%Plota ERLE
Nw = 1024;
[ERLEdB]=erle(eco, e, Nw, Fs);

%Reproduz o sinal de áudio original, o eco e o erro
soundsc(locutor,Fs);
pause(25);
soundsc(eco,Fs);
pause(25);
soundsc(e,Fs);
pause(25);