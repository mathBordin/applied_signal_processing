%%
%Experiência 5 - Codifcação de imagens segundo o padrão JPEG
%PSI3531 - Processamento de Sinais Aplicado
%Matheus Bordin Gomes

clear; close all; clc;

% Lê imagem de entrada e converte para double
gato=imread('cat.png');
gato = double(gato);

% Mostra imagem original
figure(1);
imshow(uint8(gato));
title('Imagem original');

% Mostra os componentes RGB da iamgem original
figure(2);
subplot(1,3,1); 
imshow(uint8(gato(:, :, 1)));
title('gato R');
subplot(1,3,2); 
imshow(uint8(gato(:, :, 2)));
title('gato G');
subplot(1,3,3); 
imshow(uint8(gato(:, :, 3)));
title('gato B');

% Fatores de conversão de RGB para YCbCr
alfa_r = 0.299;
alfa_g = 0.587;
alfa_b = 0.114;

% Converção da imagem de RGB para YCbCr
gato_Y = alfa_r*gato(:, :, 1) + alfa_g*gato(:, :, 2) + alfa_b*gato(:, :, 3);
gato_Cb = (1/(2*(1-alfa_b)))*(gato(:, :, 3)-gato_Y);
gato_Cr = (1/(2*(1-alfa_r)))*(gato(:, :, 1)-gato_Y);

% Mostra os componentes YCbCr da imagem original
figure(3);
subplot(1,3,1); 
imshow(uint8(gato_Y));
title('gato Y');
subplot(1,3,2); 
imshow(uint8(gato_Cb));
title('gato Cb');
subplot(1,3,3); 
imshow(uint8(gato_Cr));
title('gato Cr');

% Sub-amostra as componentes Cb e Cr por um fator 2
down_sampling_factor = 2;
Cb_sub = gato_Cb(1:down_sampling_factor:end, 1:down_sampling_factor:end);
Cr_sub = gato_Cr(1:down_sampling_factor:end, 1:down_sampling_factor:end);

% Cálculo da taxa de compressão
component_size = size(gato(:, :, 1));
original_size = 3*component_size(1)*component_size(2);
reduced_component_size = size(Cb_sub);
new_size = component_size(1)*component_size(2) + 2*reduced_component_size(1)*reduced_component_size(2);
compression = ((original_size-new_size)/original_size)*100 %#ok<NOPTS>

% Interpolação dos sinais subamostrados
h = [ 0.25 0.5 0.25 ; 0.5 1 0.5 ; 0.25 0.5 0.25];

Cb_recovered = zeros(size(gato_Y));
Cb_recovered(1:down_sampling_factor:end, 1:down_sampling_factor:end) = Cb_sub;
Cb_recovered = filter2(h, Cb_recovered);

Cr_recovered = zeros(size(gato_Y));
Cr_recovered(1:down_sampling_factor:end, 1:down_sampling_factor:end) = Cr_sub;
Cr_recovered = filter2(h, Cr_recovered);

% Mostra os componentes YCbCr recuperados
figure(4);
subplot(1,3,1); 
imshow(uint8(gato_Y));
title('gato Y');
subplot(1,3,2); 
imshow(uint8(Cb_recovered));
title('gato Cb recuperado');
subplot(1,3,3); 
imshow(uint8(Cr_recovered));
title('gato Cr recuperado');

% Conversão para componentes RGB
R_recovered = 2*(1-alfa_r)*Cr_recovered + gato_Y;
B_recovered = 2*(1-alfa_b)*Cb_recovered + gato_Y;
G_recovered = (gato_Y-alfa_r*R_recovered - alfa_b*B_recovered)/alfa_g;

% Mostra os componentes RGB da imagem recuperada
figure(5);
subplot(1,3,1); 
imshow(uint8(R_recovered));
title('gato recuperado R');
subplot(1,3,2); 
imshow(uint8(B_recovered));
title('gato recuperado G');
subplot(1,3,3); 
imshow(uint8(G_recovered));
title('gato recuperado B');

% Recupera imagem
gato_recovered = cat(3, R_recovered, G_recovered, B_recovered);

% Mostra imagem recuperada
figure(6);
imshow(uint8(gato_recovered));
title('Imagem recuperada');

% Cálculo do erro quadrático médio
MSE = 0;
for l = 1:component_size(1)
    for c = 1:component_size(2)
        MSE = MSE + (norm(gato(l,c)-gato_recovered(l,c))^2);
    end
end
MSE = MSE/(component_size(1)*component_size(2)) %#ok<NOPTS>

% Cálculo do peak signal-to-noise ratio
B=8;
PSNR = 10*log10((2^B-1)^2/MSE) %#ok<NOPTS>