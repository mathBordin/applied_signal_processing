function xq=midtreadQ(x,Nbits,fator)
% xq=midtreadQ(x,Nbits,fator)
% Essa fun��o faz a quantiza��o uniforme de x com Nbits em [-fator, +fator]

alpha=2^(Nbits-1)/fator;
xq=(floor(alpha*x+0.5))/alpha;

 
