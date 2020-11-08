%%
%Experiência 3 - Como o som é proccessado em um MP3 player
%PSI3531 - Processamento de Sinais Aplicado
%Experiência 3
%Matheus Bordin Gomes

clear; close all; clc;

%%
%Banco de filtros QMF com dois canais

%%
%Gera os filtros QMF
h0 = importdata('h0_QMF.mat'); %Lê resposta impulsiva h0
n = 1:1:length(h0);
Fa = 8000;

%Gera os filtros h1, g0 e g1
h1 = real(exp(1i*pi*n).*h0);
g0 = 2*h0;
g1 = -2*h1;
h = [h0 ; h1];
g = [g0 ; g1];

%Plota a resposta em frequência dos filtros
figure(1);
freqz(h0,1,1024,Fa);
title('Resposta em frequência de h0');
figure(2);
freqz(h1,1,1024,Fa);
title('Resposta em frequência de h1');
figure(3);
freqz(g0,1,1024,Fa);
title('Resposta em frequência de g0');
figure(4);
freqz(g1,1,1024,Fa);
title('Resposta em frequência de g1');

%%
%Gera sinal chirp linear amostrado
Fa = 8000;                              %Frequência de amostragem
Ta = 1/Fa;                              %Período de amostragem
ns = 1:1:4*Fa;
f = 2*pi*linspace(0,4000,length(ns));    
A0 = f./(2*Ta*ns);
x = cos(A0.*(Ta^2).*(ns.^2));            %Sinal chirp linear
%Gera espectrograma do sinal x (chirp linear)
figure(5);
spectrogram(x,hamming(512),256,1024,Fa,'yaxis');
title('Linear Chirp');
sound(x,Fa);
pause(4);

%%
%Filtra o sinal chirp gerado pelo banco de filtros
M = 2;
L = 2;
x_filtered = zeros(M,length(x));
v_sub = zeros(M,round(length(x)/M));
v_sup = zeros(M,length(x));
y = zeros(M,length(x));

%Primeira filtragem
for i=1:M
    x_filtered(i,:) = filter(h(i,:),1,x);
    v_sub(i,:) = x_filtered(i,1:M:end);
    v_sup(i,1:L:end) = v_sub(i,:);
    y(i,:) = filter(g(i,:),1,v_sup(i,:));
end

%Sinal recuperado
y_out = sum(y,1); 

%Espectrograma do sinal recuperado
figure(6);
spectrogram(y_out,hamming(512),256,1024,Fa,'yaxis');
title('Linear Chirp Recovered');

%Toca o sinal recuperado
sound(y_out,Fa);
pause(4);

%%
%Plota o erro de recuperacao
ordem = length(h(1,:))-1;
figure(7);
error = (x(1:(end-ordem)) - y_out((ordem+1):end));
plot(ns(1:(end-ordem))*Ta,error);
title('Erro de reconstrucao');