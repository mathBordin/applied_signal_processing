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

% Faz com que o tamanho das matrizes seja múltiplo de 6
component_size = size(gato(:, :, 1));
if mod(component_size(1),8) ~= 0
    px = 8-mod(component_size(1),8);
else
    px = 0;
end
if mod(component_size(2),8) ~= 0
    py = 8-mod(component_size(2),8);
else
    py = 0;
end
gato_Y = padarray(gato_Y, [px py], 'symmetric', 'post');
gato_Cb = padarray(gato_Cb, [px py], 'symmetric', 'post');
gato_Cr = padarray(gato_Cr, [px py], 'symmetric', 'post');

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

% Aplica DCT em blocos 8x8
dct = @(block_struct) dct2(block_struct.data,[8 8]);

gato_Y_dct = blockproc(gato_Y-128, [8 8], dct);
gato_Cb_dct = blockproc(Cb_sub-128, [8 8], dct);
gato_Cr_dct = blockproc(Cr_sub-128, [8 8], dct);

% Correlação dos coeficientes DC na componente Y (dct8x8(1,1))
figure(4);
plot(gato_Y_dct(9:8:end,1),gato_Y_dct(1:8:end-8,1),'x');
title('Correlação do nível DC em gato Y');
xlabel('GatoY_DC(n+1)');
ylabel('GatoY_DC(n)');

% Correlação dos coeficientes AC na componente Y (dct8x8(3,4))
figure(5);
plot(gato_Y_dct(3,12:8:end),gato_Y_dct(3,4:8:end-8),'x');
title('Correlação do nível AC em gato Y');
xlabel('GatoY_AC(n+1)');
ylabel('GatoY_AC(n)');

% Quantização dos coeficientes da DCT
k = 1;
Q = [ 8  11  10  16  24   40  51  61 % Tabela de quantizaçaõ JPEG
     12  12  14  19  26   58  60  55        
     14  13  16  24  40   57  69  56 
     14  17  22  29  51   87  80  62
     18  22  37  56  68  109 103  77
     24  35  55  64  81  104 113  92 
     49  64  78  87  103 121 120 101 
     72  92  95  98  112 100 103  99] * k;
quant = @(block_struct) round(block_struct.data ./ Q);

gato_Y_dct = blockproc(gato_Y_dct, [8 8], quant);
gato_Cb_dct = blockproc(gato_Cb_dct, [8 8], quant);
gato_Cr_dct = blockproc(gato_Cr_dct, [8 8], quant);

% Codificação DC
amp_dc = 2047;
prob_dc = zeros(2*amp_dc+1,1);
for i=(-1*amp_dc):1:amp_dc
    prob_dc(i+amp_dc+1)=sum(sum(gato_Y_dct(1:8:end,1:8:end) == i));
    prob_dc(i+amp_dc+1)=prob_dc(i+amp_dc+1)+sum(sum(gato_Cb_dct(1:8:end,1:8:end) == i));
    prob_dc(i+amp_dc+1)=prob_dc(i+amp_dc+1)+sum(sum(gato_Cr_dct(1:8:end,1:8:end) == i));
end
dc_size = (numel(gato_Y_dct(1:8:end,1:8:end))+2*numel(gato_Cb_dct(1:8:end,1:8:end)));
prob_dc = prob_dc/dc_size;

entropia_dc = -1*sum(prob_dc(prob_dc~=0).*log2(prob_dc(prob_dc~=0))) %#ok<NOPTS>

% Codificação AC
amp_ac = 1023;
prob_ac = zeros(2*amp_ac+1,1);
for i=(-1*amp_ac):1:amp_ac
    prob_ac(i+amp_ac+1)=sum(sum(gato_Y_dct == i))-sum(sum(gato_Y_dct(1:8:end,1:8:end) == i));
    prob_ac(i+amp_ac+1)=prob_ac(i+amp_ac+1)+sum(sum(gato_Cb_dct == i))-sum(sum(gato_Cb_dct(1:8:end,1:8:end) == i));
    prob_ac(i+amp_ac+1)=prob_ac(i+amp_ac+1)+sum(sum(gato_Cr_dct == i))-sum(sum(gato_Cr_dct(1:8:end,1:8:end) == i));
end
ac_size = ((numel(gato_Y_dct))+2*numel(gato_Cb_dct)-dc_size);
prob_ac = prob_ac/ac_size;

entropia_ac = -1*sum(prob_ac(prob_ac~=0).*log2(prob_ac(prob_ac~=0))) %#ok<NOPTS>

% Tamanho esperado do arquivo JPEG
gato_jpeg_size = entropia_dc*dc_size+entropia_ac*ac_size %#ok<NOPTS>

% Recupera componentes
rec = @(block_struct) idct2(round(block_struct.data .* Q), [8 8]);

gato_Y_idct = blockproc(gato_Y_dct, [8 8], rec) + 128;
gato_Cb_idct = blockproc(gato_Cb_dct, [8 8], rec) + 128;
gato_Cr_idct = blockproc(gato_Cr_dct, [8 8], rec) + 128;

% Interpolação dos sinais recuperados
h = [ 0.25 0.5 0.25 ; 0.5 1 0.5 ; 0.25 0.5 0.25];

Cb_rec_aux = zeros(size(gato_Y_idct));
Cb_rec_aux(1:down_sampling_factor:end, 1:down_sampling_factor:end) = gato_Cb_idct;
Cb_rec_aux = filter2(h, Cb_rec_aux);

Cr_rec_aux = zeros(size(gato_Y_idct));
Cr_rec_aux(1:down_sampling_factor:end, 1:down_sampling_factor:end) = gato_Cr_idct;
Cr_rec_aux = filter2(h, Cr_rec_aux);

Cb_recovered = Cb_rec_aux(1:end-px,1:end-py);
Cr_recovered = Cr_rec_aux(1:end-px,1:end-py);
Cy_recovered = gato_Y_idct(1:end-px,1:end-py);

% Mostra os componentes YCbCr recuperados
figure(6);
subplot(1,3,1); 
imshow(uint8(Cy_recovered));
title('gato Y');
subplot(1,3,2); 
imshow(uint8(Cb_recovered));
title('gato Cb recuperado');
subplot(1,3,3); 
imshow(uint8(Cr_recovered));
title('gato Cr recuperado');

% Conversão para componentes RGB
R_recovered = 2*(1-alfa_r)*Cr_recovered + Cy_recovered;
B_recovered = 2*(1-alfa_b)*Cb_recovered + Cy_recovered;
G_recovered = (Cy_recovered-alfa_r*R_recovered - alfa_b*B_recovered)/alfa_g;

% Mostra os componentes RGB da imagem recuperada
figure(7);
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
figure(8);
imshow(uint8(gato_recovered));
title('Imagem recuperada');

% Cálculo da razão coeficientes nulos/pixels
nulls = numel(gato_Y_dct(gato_Y_dct == 0));
nulls = nulls + numel(gato_Cb_dct(gato_Cb_dct == 0));
nulls = nulls + numel(gato_Cr_dct(gato_Cr_dct == 0));
pixels = numel(gato);
pixels_nulls = pixels/nulls  %#ok<NOPTS>

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