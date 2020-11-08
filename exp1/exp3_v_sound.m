%Experi�ncia 1 - Filtros Adaptativos
%Matheus Bordin Gomes

clear; close all; clc; 

%%%
%Cancelamento de Eco Ac ?ustico
%%%

[locutor,Fs] = audioread('locutor.wav'); %Carrega sinal de �udio
load('ri.mat'); %Carrega resposta ao impulso do ambiente que gera eco

%Gera eco do sinal de �udio com interfer�ncia
[interf, Fs_2] = audioread('eng.wav');
eco = filter(ri256, 1, locutor) + interf;

%NLMS
mi = 0.5;
delta = 10^(-5);
M = 256;

[W,e,y]=nlms(locutor, eco, M, size(locutor,1), mi, delta);

%Plota ERLE
Nw = 1024;
[ERLEdB]=erle(eco, e, Nw, Fs);

%Reproduz o sinal de �udio original, o eco e o erro
sound(locutor,Fs);
pause(25);
sound(eco,Fs);
pause(25);
sound(e,Fs);