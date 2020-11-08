%%
%Experiência 3 - Como o som é proccessado em um MP3 player
%PSI3531 - Processamento de Sinais Aplicado
%Experiência 3
%Matheus Bordin Gomes

clear; close all; clc;

%%
%Banco de filtros com dois canais

%%
%Gera sinal chirp linear amostrado
Fa = 8000;                              %Frequência de amostragem
Ta = 1/Fa;                              %Período de amostragem
n = 1:1:4*Fa;
f = 2*pi*linspace(0,4000,length(n));    
A0 = f./(2*Ta*n);

x = cos(A0.*(Ta^2).*(n.^2));            %Sinal chirp linear

%Gera espectrograma do sinal x (chirp linear)
figure(1);
spectrogram(x,hamming(512),256,1024,Fa,'yaxis');
title('Linear Chirp');

%Toca o sinal x (chirp linear)
%sound(x,Fa);
%pause(4);

%%
%Subamostragem/Superamostragem
M = 2;                  
x_sub = x(1:M:end);                     %Sinal x subamostrado

L = 2;
x_sup = zeros(1,length(n));
x_sup(1:L:end) = 2*x_sub(:);        %Sinal x_sub superamostrado com potência corrigida

%Gera o espectrograma do sinal x_sup
figure(2);
spectrogram(x_sup,hamming(512),256,1024,Fa,'yaxis');
title('Linear Chirp Supersampled');

%Toca o sinal x_sup
%sound(x_sup,Fa);
%pause(4);

%%
%Projeto do filtro passa-baixas FIR
dev = [2*10^(-4) 10^(-4)];              %Ripple na banda de passagem e de rejeição em dB
f_filter = [1960 2040];                 %Frequências de interesse
a = [1 0];                              %Amplitudes desejadas

%Projeta o filtro
[n_filter,fo,ao,w] = firpmord(f_filter,a,dev,Fa);
b_lp = firpm(n_filter,fo,ao,w);

%Plota a responsta em frequência do filtro passa-baiaxas projetado
figure(3);
freqz(b_lp,1, 1024, Fa);

%%
%Filtragem do sinal superamostrado
chirp_sup_filtered = filter(b_lp,1,x_sup); %Sinal superamostrado filtrad

%Plota o sinal chirp_sup_filtered
figure(4);
plot(n*Ta, chirp_sup_filtered);
title('Linear Chirp Supersampled after Filtering');

%Gera o espectrograma do sinal chirp_sup_filtered
figure(5);
spectrogram(chirp_sup_filtered,hamming(512),256,1024,Fa,'yaxis');
title('Linear Chirp Supersampled after Filtered');

%Toca o sinal x_sup
%sound(chirp_sup_filtered,Fa);
%pause(4);

%%
%Correção do aliasing
x_anti_aliasing = filter(b_lp,1,x);                                 %Filtra sinal original para eliniar componentes maiores do que 2 kHz
x_anti_aliasing_sub = x_anti_aliasing(1:M:end);                     %Subamostra o sinal com anti-aliasing

x_anti_aliasing_sup = zeros(1,length(n));
x_anti_aliasing_sup(1:L:end) = x_anti_aliasing_sub(:);            %Sinal superamostrado com anti-aliasing com potência corrigida

x_anti_aliasing_sup = filter(2*b_lp,1,x_anti_aliasing_sup);           %Retira componentes de alta frequência do sinal x_anti_aliasing_sup

%Gera o espectrograma do sinal chirp_sup_filtered
figure(6);
spectrogram(x_anti_aliasing_sup,hamming(512),256,1024,Fa,'yaxis');
title('Linear Chirp Filtered Supersampled - Anti Aliasing');

%Toca o sinal x_anti_aliasing_sup_filtered
%sound(x_anti_aliasing_sup,Fa);
%pause(4);

%%
%Solução para preservar altas frequências
n_filter = 1:1:length(b_lp);
b_hp = real(exp(1i*pi*n_filter).*b_lp);                               %Gera filtro passa-altas a partir do passa-baixas

%Plota a responsta em frequência do filtro passa-altas projetado
figure(7);
freqz(b_hp,1, 1024, Fa);

x_anti_aliasing_hi = filter(b_hp,1,x);                                %Elimina frequências abaixo de 2kHz do sinal original  
 
x_anti_aliasing_hi_sub = x_anti_aliasing_hi(1:M:end);                 %Sinal x_anti_aliasing_hi subamostrado

x_anti_aliasing_hi_sup = zeros(1,length(n));
x_anti_aliasing_hi_sup(1:L:end) = x_anti_aliasing_hi_sub(:);        %Sinal x_anti_aliasing_hi_sup superamostrado com potência corrigida

b_hp = 2*b_hp;
x_anti_aliasing_hi_sup = filter(b_hp,1,x_anti_aliasing_hi_sup);       %Retira componentes de alta frequência do sinal x_anti_aliasing_hi_sup

%Gera o espectrograma do sinal x_anti_aliasing_hi_sup
figure(8);
spectrogram(x_anti_aliasing_hi_sup,hamming(512),256,1024,Fa,'yaxis');
title('Linear Chirp Filtered Supersampled - Anti Aliasing (High Frequencies)');

%Toca o sinal x_anti_aliasing_hi_sup
%sound(x_anti_aliasing_hi_sup,Fa);
%pause(4);

%%
%Restaura o sinal original
x_recover = x_anti_aliasing_sup + x_anti_aliasing_hi_sup;

%Gera o espectrograma do sinal x_recover
figure(9);
spectrogram(x_recover,hamming(512),256,1024,Fa,'yaxis');
title('Linear Chirp - Recovered Signal');

%Toca o sinal x_anti_aliasing_hi_sup
%sound(x_recover,Fa);
%pause(4)

%Plota o erro de recuperacao
figure(10);
ordem = length(b_lp)-1;
error = (x(1:(end-ordem)) - x_recover(ordem+1:end));
plot(n(1:(end-ordem))*Ta,error);
title('Erro de reconstrucao');
