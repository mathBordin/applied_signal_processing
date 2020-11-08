%%
% Prova Final - PSI3531 (Processamento de Sinais Aplicado)
% Questão 2c. 
% (Como o som é proccessado em um MP3 player)
% Matheus Bordin Gomes

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
%Lê sinal de áudio usando como meio, gera sequência aleatório, calcula 
%ganho, gera filtro fir passa-baixas e gera marca d'água.
[x,Fs] = audioread('BeethovenExp4.wav'); %Sinal de áudio
R = 100;
Nb = Fs/R;

%(Questão 2c.1)
%Projeto do filtro FIR passa-baixas 
f_filtro = [9000 11000]; %Frequências em que as amplitudes são especificadas
f_filtro = f_filtro/(Fs/2); %Normaliza o vetor de frequências
f_filtro = [0 f_filtro 1];
n_filtro = 50; %Ordem do filtro
desired_amp = [1 1 0 0]; %Amplitudes especificadas para cada uma das frequências
%Retorna os coeficientes do filtro
b_fir = firpm(n_filtro,f_filtro,desired_amp);
%Plota a resposta em frequência do filtro
figure(2);
freqz(b_fir,1,1024,Fs)

%(Questão 2c.2)
%Sequência aleatória 
rand('seed',12345); 
c = 2*round(rand(Nb+(n_filtro/2),1))-1;
c = filter(b_fir,1,c);
c = c(n_filtro/2+1:end);
c = c*(sqrt(Nb)/norm(c));
var_c = var(c);

%Gera sinal modulado
v = kron(bm,c);
w = zeros(size(v));

% Calcula o sinal de marca d'água por modelagem psicoacustica
Nfiltro = 50;
Nframe = 512;
p = Nfiltro;

correlacao_v= xcorr(v,Nfiltro,'biased');

rest = mod(length(x),Nframe);
if rest == 0
    num_frames = length(x)/Nframe;
else 
    num_frames = ceil(length(x)/Nframe);
end 

mask = zeros(Nframe/2, num_frames);
b0 = zeros(num_frames, 1);
ai = zeros(Nfiltro+1, num_frames);
zi = zeros(1,p);
h_wiener=zeros(num_frames,2*Nfiltro+1);

%Filtragem com modelo pseudo-acustica
for m = 0:num_frames-2
    mask(:, m+1) = psychoacoustical_model(x(m*Nframe+1:m*Nframe+Nframe));
    mask(end/2:end, m+1) = -100; %(Questão 2c.3)
    [b0(m+1), ai(:,m+1)] = shaping_filter_design(mask(:, m+1), Nfiltro);
    [w(m*Nframe+1:m*Nframe+Nframe), zi] = filter(b0(m+1), ai(:,m+1), v(m*Nframe+1:m*Nframe+Nframe), zi);
end
mask(:,num_frames) = psychoacoustical_model(x((num_frames-1)*Nframe+1:end));
mask(end/2:end, num_frames) = -100; %(Questão 2c.3)
[b0(num_frames), ai(:,num_frames)] = shaping_filter_design(mask(:,num_frames), Nfiltro);
w((num_frames-1)*Nframe+1:end) = filter(b0(num_frames), ai(:,num_frames), v((num_frames-1)*Nframe+1:end), zi);

%Plota a PSD de um trecho do sinal de áudio
figure(3);
pwelch(x(num_frames*Nframe/2+1:num_frames*Nframe/2+Nframe),[ ],[ ],[ ],2);
title('Densidade espectral de potência do trecho numframes/2 do sinal de audio');

%Plota resposta em frequência do filtro de modelagem perceptiva do trecho
figure(4);
freqz(b0(num_frames/2+1), ai(:,num_frames/2+1));
title('Resposta em frequência do filtro de modelagem perceptiva do trecho numframes/2');

%Plota o limiar de mascaramento do trecho
figure(5);
stairs(linspace(0,1,length(mask(:, num_frames/2+1))),mask(:, num_frames/2+1));
title('Resposta em frequência do limiar de mascaramento do trecho');
ylabel('Limiar de mascaramento (dB)');
xlabel('Frequência normalizada (x pi rad/amostra)');
    
%Plota sinal de marca d'agua
figure(6);
plot((1:length(w))',w);
title('Sinal de marca dágua (sinal modulado)');
ylabel('Amplitude');
xlabel('Amostras');

%Plota sinal de audio
figure(7);
plot((1:length(x))',x);
title('Sinal de audio');
ylabel('Amplitude');
xlabel('Amostras');

%Trecho do sinal de marca d'água no tempo
figure(8);
plot((253*Nb:255*Nb+1)',w(253*Nb:255*Nb+1));
title('Sinal de marca dágua (sinal modulado) - zoom');
ylabel('Amplitude');
xlabel('Amostras');

%Densidade espectral de potência do sinal modulado
figure(9);
pwelch(w,[ ],[ ],[ ],2);
title('Densidade espectral de potência do sinal de marca dágua (sinal modulado)');

%Sinal emitido
sinal_emitido = kron(bm,ones(Nb,1));

%Trecho do sinal de emitido no tempo
figure(10);
plot((253*Nb:255*Nb+1)',sinal_emitido(253*Nb:255*Nb+1));
title('Trecho do sinal de emitido');
ylabel('Amplitude');
xlabel('Amostras');

%Densidade espectral de potência do sinal emitido
figure(11);
pwelch(sinal_emitido,[ ],[ ],[ ],2);
title('Densidade espectral de potência do sinal emitido');

%%
%Gera sinal de saída
%(Questão 2c.4)
y = x + w;
yc = mp3_codec(y,Fs);
yc = filter(b_fir,1, [yc ; zeros(n_filtro/2,1)]);
yc = yc(n_filtro/2+1:end);

%Plota sinais
figure(12);
plot((253*Nb:255*Nb+1)',w(253*Nb:255*Nb+1), 'g',(253*Nb:255*Nb+1)',x(253*Nb:255*Nb+1), 'b',(253*Nb:255*Nb+1)',yc(253*Nb:255*Nb+1),'r');
title('Sinais de marca dágua, de áudio e da soma');
ylabel('Amplitude');
xlabel('Amostras');
legend('Sinal de marca dágua', 'Sinal de áudio', 'Sinal somado');

%%
%Receptor

%%
% Recupera sinal de marca d'agua com o inverso do filtro com a estimação da
% mascara da modelagem pseudo-acustica
mask_y = zeros(Nframe/2, num_frames);
b0_y = zeros(num_frames, 1);
ai_y = zeros(Nfiltro+1, num_frames);
zi_y = zeros(1,p);
z = zeros(size(yc));

for m = 0:num_frames-2
    mask_y(:, m+1) = psychoacoustical_model(yc(m*Nframe+1:m*Nframe+Nframe));
    mask_y(end/2:end, m+1) = -100; %(Questão 2c.4)
    [b0_y(m+1), ai_y(:,m+1)] = shaping_filter_design(mask_y(:, m+1), Nfiltro);
    [yc(m*Nframe+1:m*Nframe+Nframe), zi_y] = filter(ai_y(:,m+1), b0_y(m+1), yc(m*Nframe+1:m*Nframe+Nframe), zi_y);
    
    corr= xcorr(yc(m*Nframe+1:m*Nframe+Nframe),2*Nfiltro,'biased');
    h_wiener(m+1,:) = linsolve(toeplitz(corr(101:end)),correlacao_v);
    z(m*Nframe+1:m*Nframe+Nframe) = filter(h_wiener(m+1,:),1,yc(m*Nframe+1:m*Nframe+Nframe));
    pot_z = var(z(m*Nframe+1:m*Nframe+Nframe));
    if pot_z ~= 0
        z(m*Nframe+1:m*Nframe+Nframe) = z(m*Nframe+1:m*Nframe+Nframe)/sqrt(pot_z);
    end
end
mask_y(:,num_frames) = psychoacoustical_model(yc((num_frames-1)*Nframe+1:end));
mask_y(end/2:end, num_frames) = -100; %(Questão 2c.4)
[b0_y(num_frames), ai_y(:,num_frames)] = shaping_filter_design(mask_y(:,num_frames), Nfiltro);    
yc((num_frames-1)*Nframe+1:end) = filter(ai_y(:,num_frames),b0_y(m+1), yc((num_frames-1)*Nframe+1:end), zi_y);

corr= xcorr(yc((num_frames-1)*Nframe+1:end),2*Nfiltro,'biased');
h_wiener(num_frames,:) = linsolve(toeplitz(corr(101:end)),correlacao_v);
z((num_frames-1)*Nframe+1:end) = filter(h_wiener(num_frames,:),1,yc((num_frames-1)*Nframe+1:end));
pot_z = var(z(m*Nframe+1:m*Nframe+Nframe));
if pot_z ~= 0
    z(m*Nframe+1:m*Nframe+Nframe) = z(m*Nframe+1:m*Nframe+Nframe)/sqrt(pot_z);
end
z = [z(Nfiltro+1:end) ; zeros(Nfiltro,1)];

%%
%Calcula produto escalar normalizado entre o sinal recebido e a sequência
%pseudo-aletória e decide bits recebidos
alfa_m = zeros(M,1);
beta_m = zeros(M,1);
for i = 1:M
    %Produto escalar normalizado
    alfa_m(i) = z((i-1)*Nb+1:(i-1)*Nb+Nb)'*c/Nb;
    %Decisão de bits recebidos
    if alfa_m(i)<=0
        beta_m(i) = -1;
    else
        beta_m(i) = 1;
    end
end

%Plota o gráfico dos bits emitidos, dos bits recuperados e do produto
%escalar normalizado
figure(13);
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
figure(14);
imshow(img_rcv,[ ]);
title('Imagem Recuperada');

%Calcula BER
BER = sum(abs(str2num(bin_signal)-bits_rcv))/length(bin_signal) %#ok<ST2NM,NOPTS>

%Calcula índice de similaridade estrutural
mssim = ssim_index(img, img_rcv,[0.01 0.03], fspecial('gaussian', 11, 1.5), 255) %#ok<NOPTS>

% sound(yc,Fs); %Descomente para ouvir o sinal com a marca d'água