function [v1]=cheb_cumModeCell_rotating(t,m,wm,OmagR,r,Gw1,Gw2,Ratio,lx1,lx2)
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
        [Gww]=cumModeCell(m,Gww,Gw1,kkk,t,solutime,wm,OmagR);   
        v1{solutime}(:,:,kkk)=real(Gww(r));
     end
end
toc         
         
         
%          tsignal.cubes(solutime).v(1,:,:,kkk)=real(Gww{solutime}{1,kkk}(r));
%          tsignal.cubes(solutime).v(2,:,:,kkk)=real(Tmm1{solutime}{1,kkk}(r));
%          tsignal.cubes(solutime).v(3,:,:,kkk)=real(Tmm2{solutime}{1,kkk}(r));
%          tsignal.cubes(solutime).v(4,:,:,kkk)=real(Tmm3{solutime}{1,kkk}(r));
   
%          GGw{solutime}(:,:,kkk)=Gww{1,kkk}(r);
%          TTm1{solutime}(:,:,kkk)=Tmm1{1,kkk}(r);
%          TTm2{solutime}(:,:,kkk)=Tmm2{1,kkk}(r);
%          TTm3{solutime}(:,:,kkk)=Tmm3{1,kkk}(r);

