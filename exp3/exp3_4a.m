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
sound(x,Fa);
pause(2);

%%
%Filtra o sinal chirp gerado pelo banco de filtros
M = 32;
L = 32;
x_filtered = zeros(M,length(x));
v_sub = zeros(M,ceil(length(x)/M));
v_sup = zeros(M,length(x));
y = zeros(M,length(x));

n_bits = 4;
fator_q = 1;
x_in = midtreadQ(x,n_bits,fator_q);

%Primeira filtragem
for i=1:M
    x_filtered(i,:) = filter(PQMF32_Hfilters(i,:),1,x_in);
    v_sub(i,:) = x_filtered(i,1:M:end);
    v_sup(i,1:L:end) = v_sub(i,:);
    y(i,:) = filter(PQMF32_Gfilters(i,:),1,v_sup(i,:));
    figure(3);
    plot(ns*Ta,y(i,:));
    title('Saída dos filtros G');
    hold on;
end

%Sinal recuperado
y_out = sum(y,1); 

%Espectrograma do sinal recuperado
figure(4);
spectrogram(y_out,hamming(512),256,1024,Fa);
title('Violin Recovered');

%Toca o sinal recuperado
sound(y_out,Fa);
pause(2);

%%
%Calcula o erro de recuperação 
ordem = length(PQMF32_Hfilters(1,:))-1;
error = x(1:(end-ordem))- y_out(ordem+1:end);

figure(5);
plot(ns(1:(end-ordem))*Ta,error);
title('Erro de recuperação');

SNR_PQMF=snr(x(1:(end-ordem)), y_out(ordem+1:end),0)