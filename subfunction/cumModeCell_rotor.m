function [Gww]=cumModeCell_rotor(m,Gww,Gw1,kkk,t,solutime,w)        
for kk=1:length(m)
             %Pww{1,kkk}=Pww{1,kkk}+Pw{1,kk}{1,kkk};%��m����
            Gww=Gww+Gw1{1,kk}{1,kkk}*exp(-sqrt(-1)*(w(kk))*t(solutime));%��m����     
end