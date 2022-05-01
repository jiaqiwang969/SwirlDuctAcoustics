function [Gww]=cumModeCell(m,Gww,Gw1,kkk,t,solutime,wm,OmagR)        
for kk=1:length(m)
             %Pww{1,kkk}=Pww{1,kkk}+Pw{1,kk}{1,kkk};%对m积分
            Gww=Gww+Gw1{1,kk}{1,kkk}*exp(-sqrt(-1)*(wm+m(kk)*OmagR)*t(solutime));%对m积分     
end