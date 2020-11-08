function [y]=lappedQadap(x,fa,H,G,Nbits)
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
% Nbits: N�mero de bits usados na quantiza��o

x=x(:)';
Nx=length(x);
[M,Nh]=size(H);

output_signal=zeros(1,Nx);
Nframes=fix((Nx-Nh+M)/M);
Nframes=fix(Nframes/12)*12;
subbands=zeros(Nframes,M);
subbandsQ=zeros(Nframes,M);

for i=1:Nframes
    input_frame=x((i-1)*M+1:(i-1)*M+Nh);
    subbands(i,:)=(G*input_frame')';
end
for i=1:12:Nframes
    fator= max(abs(subbands(i:i+11,:)));
    for k=1:M
        subbandsQ(i:i+11,k) = midtreadQ(subbands(i:i+11,k),Nbits,fator(k));
    end
end
for i=1:Nframes    
    output_frame = G'*subbandsQ(i,:)';
    output_signal((i-1)*M+1:(i-1)*M+Nh)= output_signal((i-1)*M+1:(i-1)*M+Nh)+output_frame';
end
y=output_signal(Nh:end-Nh);
