%%
%Experiência 3 - Como o som é proccessado em um MP3 player
%PSI3531 - Processamento de Sinais Aplicado
%Experiência 3
%Matheus Bordin Gomes

clear; close all; clc;

%%
%Codificação de audio perceptual
load('Analise_Sintese32.mat'); %Lê banco de filtros
Fa = 44100;

%Plota resposta em frequência dos filtos h
figure(1);
for i=1:32
    [H,W] = freqz(PQMF32_Hfilters(i,:),1,1024,Fa);
    plot(W, 20*log10(abs(H)));
    hold on;
end
title('Resposta em frequência do banco PQMF32-H');

%%
%Gera sinal chirp linear (0 a 22050Hz)
Ta = 1/Fa;                              %Período de amostragem
[x_t,Fa] = audioread('violino.wav');
x = x_t';
ns = 1:1:length(x);

%Gera espectrograma do sinal x (chirp linear)
figure(2);
spectrogram(x,hamming(512),256,1024,Fa);
title('Violin');

%Toca o sinal x
%sound(x,Fa);
%pause(2);

%Plota o sinal x
figure(6);
plot(x);
title('original');

%%
%Filtra o sinal chirp gerado pelo banco de filtros
n_bits = 4;
y_out=lappedQadap(x, Fa, PQMF32_Hfilters, PQMF32_Gfilters, n_bits);

%Espectrograma do sinal recuperado
figure(4);
spectrogram(y_out,hamming(512),256,1024,Fa);
title('Violin Recovered');

%Toca o sinal recuperado
%sound(y_out,Fa);
%pause(2);

%%
%Calcula o erro de recuperação 
ordem = length(PQMF32_Hfilters(1,:))-1;
error = x(ordem+1:end-ordem-1)- y_out(1:end);

figure(5);
plot(ns(ordem+1:end-ordem-1)*Ta,error);
title('Erro de recuperação');

SNR_PQMF=snr(x(ordem+1:end-ordem-1), y_out(1:end),0)