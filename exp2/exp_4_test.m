%Experiência 2 - Codificação de sinais de voz por análise preditiva
%Matheus Bordin Gomes

clear; close all; 

%Leitura do sinal de áudio
[sinal, fs] = audioread('antarctica.wav');

%Recorte do sinal de áudio
trecho = sinal(200:439);

%Cálculo dos parâmetros
[ak, sig2]=lpc(trecho.*hamming(240), 10);
freqz(sqrt(sig2), ak, 512); %Como se calcula o ganho?
hold on;
periodogram(trecho,[],512);

%Estimador do pitch
%pitch=yaapt(sinal,fs,1,[],1,1);
%pitch = [0 pitch];