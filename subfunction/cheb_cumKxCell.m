function [GNk,TGm1,TGm2,TGm3]=cheb_cumKxCell(G_nm,Tgm1,Tgm2,Tgm3,Ratio,lx,lm)
%initialization
%PNk=[];
GNk=[];TGm1=[];TGm2=[];TGm3=[];
for kkk=1:lx
    %PNk{1,kkk}=chebfun(0,[Ratio, 1]);
    GNk{1,kkk}=chebfun(0,[Ratio, 1]);
    TGm1{1,kkk}=chebfun(0,[Ratio, 1]);
    TGm2{1,kkk}=chebfun(0,[Ratio, 1]);
    TGm3{1,kkk}=chebfun(0,[Ratio, 1]);
    
    for kk=1:lm
         %PNk{1,kkk}=PNk{1,kkk}+p_nm{kk,kkk};%k-Integral
         GNk{1,kkk}=GNk{1,kkk}+G_nm{kk,kkk};%k-Integral
         TGm1{1,kkk}=TGm1{1,kkk}+Tgm1{kk,kkk};%k-Integral
         TGm2{1,kkk}=TGm2{1,kkk}+Tgm2{kk,kkk};%k-Integral
         TGm3{1,kkk}=TGm3{1,kkk}+Tgm3{kk,kkk};%k-Integral
    end
end
