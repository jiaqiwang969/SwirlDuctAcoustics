function [GGw,TTm1,TTm2,TTm3]=cheb_cumModeCell(r,Gw,Tm1,Tm2,Tm3,Ratio,lx)
% tic
%  for kk=1:length(Gw)
%      for kkk=1:lx
%      Pww{1,kk}{1,kkk}=Pw{1,kk}{1,kkk}(r);    
%      Gww{1,kk}{1,kkk}=Gw{1,kk}{1,kkk}(r);    
%      Tmm{1,kk}{1,kkk}=Tm{1,kk}{1,kkk}(r);    
%      end
%  end
%  PPw=reshape(sum(cell2mat(cat(3,Pww{:})),3),length(r),360,lx);
%  GGw=reshape(sum(cell2mat(cat(3,Gww{:})),3),length(r),360,lx);
%  TTm=reshape(sum(cell2mat(cat(3,Tmm{:})),3),length(r),360,lx);
% toc
tic 
%PPw=[];
GGw=[];TTm1=[];TTm2=[];TTm3=[];
for kkk=1:lx
    %Pww{1,kkk}=chebfun(0,[Ratio, 1]);
    Gww{1,kkk}=chebfun(0,[Ratio, 1]);
    Tmm1{1,kkk}=chebfun(0,[Ratio, 1]);
    Tmm2{1,kkk}=chebfun(0,[Ratio, 1]);
    Tmm3{1,kkk}=chebfun(0,[Ratio, 1]);
    for kk=1:length(Gw)
         %Pww{1,kkk}=Pww{1,kkk}+Pw{1,kk}{1,kkk};%k-Integral
         Gww{1,kkk}=Gww{1,kkk}+Gw{1,kk}{1,kkk};%k-Integral
         Tmm1{1,kkk}=Tmm1{1,kkk}+Tm1{1,kk}{1,kkk};%k-Integral
         Tmm2{1,kkk}=Tmm2{1,kkk}+Tm2{1,kk}{1,kkk};%k-Integral
         Tmm3{1,kkk}=Tmm3{1,kkk}+Tm3{1,kk}{1,kkk};%k-Integral
    end
         %PPw(:,:,kkk)=Pww{1,kkk}(r);
         GGw(:,:,kkk)=Gww{1,kkk}(r);
         TTm1(:,:,kkk)=Tmm1{1,kkk}(r);
         TTm2(:,:,kkk)=Tmm2{1,kkk}(r);
         TTm3(:,:,kkk)=Tmm3{1,kkk}(r);
end
toc