%Experiência 2 - Codificação de sinais de voz por análise preditiva
%Matheus Bordin Gomes

clear; close all; clc;

%Leitura do sinal de áudio
[sinal, fs] = audioread('antarctica.wav');

%Deixa o sinal alinhado em 240 amostras e faz o recorte
%sinal = cat(1,sinal, zeros(240 - mod(length(sinal),240),1));
n_trechos = ceil(length(sinal)/240);
trechos = zeros(240, n_trechos);

old_index = 1;
next_index = 240;
for i = 1:n_trechos
    if(i~=n_trechos)
        trechos(:,i) = sinal(old_index:next_index);
        old_index = next_index + 1;
        next_index = old_index + 239;
    else
        trechos(1:(length(sinal)-old_index+1),i) = sinal(old_index:end);
    end
end

%Encontra pitch para trechos de 30 ms
pitch = yaapt(sinal,fs,1,[],0,1);
pitch = [0 pitch];

%Encontra os parâmetros do lpc
ak = zeros(11, n_trechos);
sig2 = zeros(n_trechos,1);

for i = 1:n_trechos
    [ak(:,i), sig2(i)] = lpc(trechos(:,i).*hamming(240), 10);
    for j = 1:11
        ak(j,i) = quantize3(ak(j,i),7);
    end
end

%Reconstrói o sinal
noise_exc = randn(80,1);
lpc_out = [];
periodos = round(fs./pitch);
for i = 1:length(periodos)
    periodos(i) = quantize3(periodos(i),5);
end

for li = 1:length(pitch)
    if pitch(li) > 0
        ones_exc = ones(80,1);
        son_exc = upsample(ones_exc,periodos(1,li)); 
        son_exc = son_exc(1:80);
        
        %gain = sqrt(periodos(1,li)*sig2(round(li/3,0)+1));
        gain = quantize3(sqrt(periodos(1,li)*sig2(round(li/3,0)+1)), 5);
        
        y_out = filter(gain, ak(:,round(li/3,0)+1),son_exc);
    else
        %gain = sqrt(sig2(round(li/3,0)+1));          
        gain = quantize3(sqrt(sig2(round(li/3,0)+1)),5); 
        
        y_out = filter(gain, ak(:,round(li/3,0)+1),noise_exc);   
    end
    
    lpc_out = horzcat(lpc_out, y_out');
end

lpc_out(isnan(lpc_out))=0;

sound(lpc_out, fs);