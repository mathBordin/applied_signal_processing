function [y]=lappedPsico(x,fa,H,G,bitrate)
% Bancos de filtros e transformadas Lapped
%
% y:  sinal de sa�da
%
% x:  sinal de entrada
% fa: frequ�ncia de amostragem
% H:  matriz M x Nh dos filtros de an�lise. Cada linha deve conter a resposta
%     impulsiva dos filros de an�lise. Assim, na linha 1, teremos a resposta
%     impulsiva do filtro da Banda 0. Na linha M, a resposta impulsiva do
%     filtro da Banda M-1. Nh � o comprimento da reposta impulsiva de cada
%     filtro
% G:  matriz M x Nh dos filtros de s�ntese. Segue o mesmo padr�o da matrix
%     H.
% bitrate: taxa de bits desej�vel

x=x(:);
Nx=length(x);
[M,Nh]=size(H);

output_signal=zeros(1,Nx);
Nframes=fix((Nx-Nh+M)/M);
Nframes=fix(Nframes/12)*12;
subbands=zeros(Nframes,M);
subbands_out=zeros(Nframes,M);
subbandsQ=zeros(Nframes,M);

for i=1:Nframes
    input_frame=x((i-1)*M+1:(i-1)*M+Nh);
    subbands(i,:)=(G*input_frame)';
end

j=1;
SMR = zeros(ceil(Nframes/12),32);
Nbits = zeros(ceil(Nframes/12),32);
SNR = zeros(ceil(Nframes/12),32);
for i=1:12:Nframes
    fator=max(abs(subbands(i:i+11,:)));
    % O frame de entrada para o modelo psicoac�stico � atrasado de 176
    % amostras, j� que ele deve estar centrado no meio do bloco local das 12
    % amostras das sub-bandas e a amostra da primeira sub-banda corresponde
    % � amostra original 256 (na verdade, corresponde � amostra 512, mas
    % os filtros de an�lise e s�ntese introduzem um atraso de 256
    % amostras). Dessa forma, o centro do primeiro frame deve estar na
    % amostra 256+11*32/2=432 e o come�o desse frame cai na amostra
    % 11*32/2=176.
    frame=x(176+(i-1)*M:176+(i-1)*M+Nh-1); 
    SMR(j,:) = MPEG1_psycho_acoustic_model1(frame);
    Nbits(j,:) = MPEG1_bit_allocation(SMR(j,:), bitrate);
    for k=1:M
        if Nbits(j,k)~=0
            subbandsQ(i:i+11,k) = midtreadQ(subbands(i:i+11,k),Nbits(j,k),fator(k));
        else
            subbandsQ(i:i+11,k) = 0;
        end
    end
    j=j+1;
end

for i=1:Nframes    
    output_frame = G'*subbandsQ(i,:)';
    output_signal((i-1)*M+1:(i-1)*M+Nh)= output_signal((i-1)*M+1:(i-1)*M+Nh)+output_frame';
end
y=output_signal(Nh:end-Nh);

figure(4);
mesh(SMR);
colormap(jet);
title('SMR(dB) x Bandas x frames');
ylabel('Bloco (frames)');
xlabel('Banda de frequ�ncia');
zlabel('SMR (dB)');

figure(5);
mesh(Nbits);
colormap(jet);
title('N�mero de bits alocados x Bandas x frames');
ylabel('Bloco (frames)');
xlabel('Banda de frequ�ncia');
zlabel('N bits');