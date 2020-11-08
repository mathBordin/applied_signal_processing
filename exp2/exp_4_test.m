%Experi�ncia 2 - Codifica��o de sinais de voz por an�lise preditiva
%Matheus Bordin Gomes

clear; close all; 

%Leitura do sinal de �udio
[sinal, fs] = audioread('antarctica.wav');

%Recorte do sinal de �udio
trecho = sinal(200:439);

%C�lculo dos par�metros
[ak, sig2]=lpc(trecho.*hamming(240), 10);
freqz(sqrt(sig2), ak, 512); %Como se calcula o ganho?
hold on;
periodogram(trecho,[],512);

%Estimador do pitch
%pitch=yaapt(sinal,fs,1,[],1,1);
%pitch = [0 pitch];