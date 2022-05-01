function [v1]=cheb_cumModeCell_rotor(t,m,w,r,Gw1,Gw2,Ratio,lx1,lx2)
for k1=1:lx2
    for k2=1:length(m)
Gw1{1,k2}{1,lx1+k1}=Gw2{1,k2}{1,k1};
    end
end
tic 
Gww=[];
parfor solutime=1:length(t)
    for kkk=1:(lx1+lx2)
        Gww=chebfun(0,[Ratio, 1]);
        [Gww]=cumModeCell_rotor(m,Gww,Gw1,kkk,t,solutime,w);   
        v1{solutime}(:,:,kkk)=real(Gww(r));
     end
end
toc         
         




