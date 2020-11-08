function [y]=lappedQ(x, fa, H, G,Nbits)
% Bancos de filtros e transformadas Lapped
%
% y:  sinal de saída
%
% x:  sinal de entrada
% fa: frequência de amostragem
% H:  matriz M x Nh dos filtros de análise. Cada linha deve conter a resposta
%     impulsiva dos filros de análise. Assim, na linha 1, teremos a resposta
%     impulsiva do filtro da Banda 0. Na linha M, a resposta impulsiva do
%     filtro da Banda M-1. Nh é o comprimento da reposta impulsiva de cada
%     filtro
% G:  matriz M x Nh dos filtros de síntese. Segue o mesmo padrão da matrix
%     H.
% Nbits: Número de bits usado na quantização


x=x(:)';
Nx=length(x);
[M,Nh]=size(H);

output_signal=zeros(1,Nx);
Nframes=fix((Nx-Nh+M)/M);
Nframes=fix(Nframes/12)*12;
subbands=zeros(Nframes,M);

for i=1:Nframes
    
    % Overlap input_frames (column vectors)
    input_frame=x((i-1)*M+1:(i-1)*M+Nh);

    % Analysis filters and downsampling
    % Since PQMF H filters are the time-reversed G filters, we use the G
    % filters matrix to simulate analysis filtering
    
    subbands(i,:) = midtreadQ((G*input_frame')',Nbits,1);
    
    % Synthesis filters
    output_frame = G'*subbands(i,:)';

    % Overlap output_frames (with delay of 511 samples)
    output_signal((i-1)*M+1:(i-1)*M+Nh)= output_signal((i-1)*M+1:(i-1)*M+Nh)+output_frame';
    
end
y=output_signal(Nh:end-Nh+1);
