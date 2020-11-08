function [ERLEdB]=erle(d,e,Nw,fa);
%
% Função para cálcudo do ERLE 
% (echo return loss enhancement)
% após o cancelador de eco
% 
% Saída 
% Estimativa do ERLE em dB
%
% Entradas
% d: sinal de eco (desejado)
% e: sinal de eco residual (erro do filtro adaptativo)
% Nw: número de amostras na janela para estimativa do eco ao longo do tempo
% fa: frequência de amostragem 
%
%

% MTMS, 08/2018

N=length(d);
Nb=floor(N/Nw);
ERLEx=zeros(1,Nb);
Ta=1/fa;

for i=1:Nb
   l=Nw*(i-1)+1:Nw*i;
   ERLEx(i)=mean(d(l).^2)/mean((e(l)+eps).^2);
end
ERLEdB=10*log10(ERLEx);

figure
t=0:(N-1)*Ta/(N/Nw-1):(N-1)*Ta;
plot(t,ERLEdB)
grid on
hold on
ylabel('     Eco                                      ERLE (dB)')
xlabel('tempo (s)')
td=t(1):(t(end)-t(1))/(length(d)-1):t(end);
plot(td,10*d/max(abs(d))-11)
plot([td(1) td(end)],[-0.5 -0.5],'k','LineWidth',1)
set(gca,'YTick',[0:10:max(ERLEdB), ceil(max(ERLEdB))]');
axis([td(1) td(end) -21 ceil(max(ERLEdB))])
