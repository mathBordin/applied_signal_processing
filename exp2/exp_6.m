%Experiência 2 - Codificação de sinais de voz com realimentação (análise
%por síntese)
%Matheus Bordin Gomes

clear; close all; clc;

%Leitura do sinal de áudio
[sinal, fs] = audioread('antarctica.wav');

%Deixa o sinal alinhado em 240 amostras e faz o recorte
sinal = cat(1,sinal, zeros(240 - mod(length(sinal),240),1));

N = 80;
L=240;
n_trechos = ceil(length(sinal)/N);
trechos = zeros(240, n_trechos);

old_index = 1;
next_index = L;
for i = 1:n_trechos
    if(i~=n_trechos)
        trechos(:,i) = sinal(old_index:next_index);
        old_index = next_index + 1;
        next_index = old_index + N - 1;
    else
        trechos(1:(length(sinal)-old_index+1),i) = sinal(old_index:end);
    end
end

%Define a base de funções aleatórias
K = 2;
Qn = 512;
Q = randn(N,Qn);

%Variáveis para guardar as condições final/inicial do filtro de trato vocal para cada
%quadro no codificador/decodificador
p = 10;
zi = zeros(1,p);
zs = zeros(1,p);

%Codificador
ak = zeros(p+1, n_trechos);
sig2 = zeros(n_trechos,1);
ganhos = zeros(K,n_trechos);
indices = zeros(K,n_trechos);
d = zeros(N,n_trechos);

for i = 1:n_trechos
    [ak(:,i), sig2(i)] = lpc(trechos(:,i).*hamming(240), p);
    sub_quadro = trechos(N+1:2*N,i);
    Qy = zeros(N,Qn);
    for j = 1:Qn
        Qy(:,j) = filter(1,ak(:,i),Q(:,j));
    end
    [y0,zi]=filter(1,ak(:,i),zeros(N,1),zi);
    e0=sub_quadro-y0;
    [ganhos(:,i), indices(:,i)] = find_Nbest_components(e0, Qy, K);
    for b = 1:K
        d(:,i) = d(:,i) + ganhos(b,i)*Q(:,indices(b,i));
    end
end

%Decodificador
sinal_saida = [];

for li = 1:n_trechos
    [y_out, zs] = filter(1, ak(:,li),d(:,li),zs); 
    y_out(isnan(y_out))=0;
    zs(isnan(zs))=0;   
    sinal_saida = horzcat(sinal_saida, y_out');
end

sound(sinal_saida, fs);