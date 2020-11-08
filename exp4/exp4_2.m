%%
%Experiência 3 - Como o som é proccessado em um MP3 player
%PSI3531 - Processamento de Sinais Aplicado
%Experiência 4
%Matheus Bordin Gomes

clear; close all; clc;

%%
%Transmissor

%%
%Lê imagem que será utilizada como sequência de bits da marca d'água e gera
%sequência de bits mapeada para{-1;+1}.
img = imread('lenna64.jpg');

%Exibe imagem transmitida
figure(1);
imshow(img,[ ]);
title('Imagem Transmitida');

%Gera sequência de bits
bin_signal_matrix = dec2bin(img);
bin_signal = reshape(bin_signal_matrix,[],1);
bm = 2*(str2num(bin_signal))-1; %#ok<ST2NM>
M = length(bm);

%%
%Lê sinal de áudio usado como meio, gera sequência aleatório, calcula 
%ganho e gera marca d'água.
[x,Fs] = audioread('BeethovenExp4.wav'); %Sinal de áudio
R = 100;
Nb = Fs/R;

%Sequência aleatória 
%rand('seed',12345);
rng(12345);
c = (2*round(rand(Nb,1))-1);
var_c = var(c);

%Gera sinal modulado
v = kron(bm,c);

% Calcula ganho por frame do sinal modulado
Dg=0.005;
g = zeros(M, 1);
for m = 0:M-1
    beta=x(m*Nb+1:m*Nb+Nb)'*c/Nb;
    if bm(m+1)==1
        if beta>=Dg
            g(m+1)=0;
        else
            g(m+1)= Dg-beta;
        end
    else
        if beta<-Dg
            g(m+1)= 0;
        else
            g(m+1)=Dg+beta;
        end
    end
    v(m*Nb+1:m*Nb+Nb)= g(m+1)*v(m*Nb+1:m*Nb+Nb);
    % Esse é o sinal de marca d'água
end

%Plota sinal de marca d'agua
figure(2);
plot((1:length(v))',v);
title('Sinal de marca dágua (sinal modulado)');
ylabel('Amplitude');
xlabel('Amostras');

%Plota ganho por frame
figure(3);
plot((1:length(g))',g);
title('Ganho por frame');
ylabel('Amplitude');
xlabel('frames');

%Plota sinal de audio
figure(4);
plot((1:length(x))',x);
title('Sinal de audio');
ylabel('Amplitude');
xlabel('Amostras');

%Trecho do sinal de marca d'água no tempo
figure(5);
plot((253*Nb:255*Nb+1)',v(253*Nb:255*Nb+1));
title('Sinal de marca dágua (sinal modulado) - zoom');
ylabel('Amplitude');
xlabel('Amostras');

%Densidade espectral de potência do sinal modulado
figure(6);
pwelch(v,[ ],[ ],[ ],2);
title('Densidade espectral de potência do sinal de marca dágua (sinal modulado)');

%Sinal emitido
sinal_emitido = kron(bm,ones(Nb,1));

%Trecho do sinal de emitido no tempo
figure(7);
plot((253*Nb:255*Nb+1)',sinal_emitido(253*Nb:255*Nb+1));
title('Trecho do sinal de emitido');
ylabel('Amplitude');
xlabel('Amostras');

%Densidade espectral de potência do sinal emitido
figure(8);
pwelch(sinal_emitido,[ ],[ ],[ ],2);
title('Densidade espectral de potência do sinal emitido');

%%
%Gera sinal de saída
y = x + v;
figure(9);
plot((253*Nb:255*Nb+1)',v(253*Nb:255*Nb+1), 'g',(253*Nb:255*Nb+1)',x(253*Nb:255*Nb+1), 'b',(253*Nb:255*Nb+1)',y(253*Nb:255*Nb+1),'r');
title('Sinais de marca dágua, de áudio e da soma');
ylabel('Amplitude');
xlabel('Amostras');
legend('Sinal de marca dágua', 'Sinal de áudio', 'Sinal somado');

%%
%Receptor

%%
%Calcula produto escalar normalizado entre o sinal recebido e a sequência
%pseudo-aletória e decide bits recebidos
alfa_m = zeros(M,1);
beta_m = zeros(M,1);
for i = 1:M
    %Produto escalar normalizado
    alfa_m(i) = y((i-1)*Nb+1:(i-1)*Nb+Nb)'*c/Nb;
    %Decisão de bits recebidos
    if alfa_m(i)<=0
        beta_m(i) = -1;
    else
        beta_m(i) = 1;
    end
end

%Plota o gráfico dos bits emitidos, dos bits recuperados e do produto
%escalar normalizado
figure(10);
graphs = [bm(1:50) alfa_m(1:50) beta_m(1:50)];
h_stairs = stairs(graphs);
h_stairs(1).Color = 'r';
h_stairs(2).Color = 'g';
h_stairs(3).Color = 'b';
title('Comparaçao entre os bits emitidos, alfa_m e bits recuperados');
legend('Bits emitidos', 'Alfa_m', 'Bits recuperados');
xlabel('Amostras');
ylabel('Amplitude');

%%
%Recupera imagem e calcula taxa de erro
bits_rcv = (beta_m+1)/2;
img_rcv = reshape(bin2dec(num2str(reshape(bits_rcv,size(bin_signal_matrix)))),size(img));

%Exibe imagem recuperada
figure(11);
imshow(img_rcv,[ ]);
title('Imagem Recuperada');

%Calcula BER
BER = sum(abs(str2num(bin_signal)-bits_rcv))/length(bin_signal) %#ok<ST2NM,NOPTS>

%Calcula índice de similaridade estrutural
mssim = ssim_index(img, img_rcv,[0.01 0.03], fspecial('gaussian', 11, 1.5), 255) %#ok<NOPTS>