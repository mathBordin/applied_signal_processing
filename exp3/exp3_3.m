%%
%Experiência 3 - Como o som é proccessado em um MP3 player
%PSI3531 - Processamento de Sinais Aplicado
%Experiência 3
%Matheus Bordin Gomes

clear; close all; clc;

%%
%Banco de filtros pseudo-QMF com 32 canais
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
ns = 1:1:20*Fa;
f = 2*pi*linspace(0,Fa/2,length(ns));    
A0 = f./(2*Ta*ns);
x = cos(A0.*(Ta^2).*(ns.^2));            %Sinal chirp linear
%Gera espectrograma do sinal x (chirp linear)
figure(2);
spectrogram(x,hamming(512),256,1024,Fa);
title('Linear Chirp');
%Toca o sinal x
%sound(x,Fa);
%pause(20);

%%
%Filtra o sinal chirp gerado pelo banco de filtros
M = 32;
L = 32;
x_filtered = zeros(M,length(x));
v_sub = zeros(M,round(length(x)/M));
v_sup = zeros(M,length(x));
y = zeros(M,length(x));

%Primeira filtragem
for i=1:M
    x_filtered(i,:) = filter(PQMF32_Hfilters(i,:),1,x);
    v_sub(i,:) = x_filtered(i,1:M:end);
    v_sup(i,1:L:end) = v_sub(i,:);
    y(i,:) = filter(PQMF32_Gfilters(i,:),1,v_sup(i,:));
    figure(3);
    plot(y(i,:));
    title('Saída dos filtros G');
    hold on;
end

%Sinal recuperado
y_out = sum(y,1); 

%Espectrograma do sinal recuperado
figure(4);
spectrogram(y_out,hamming(512),256,1024,Fa);
title('Linear Chirp Recovered');

%Toca o sinal recuperado
%sound(y_out,Fa);
%pause(20);

%%
%Calcula o erro de recuperação 
ordem = length(PQMF32_Hfilters(1,:))-1;
error = x(1:(end-ordem))- y_out((ordem+1):end);

figure(5);
plot(ns(1:(end-ordem))*Ta,error);
title('Erro de recuperação');

SNR_PQMF=snr(x(1:(end-ordem)), y_out(ordem+1:end),0)