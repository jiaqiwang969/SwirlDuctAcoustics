function [Gw,Tm1,Tm2,Tm3,Gwn,TGmn1,TGmn2,TGmn3]=greenfun_span2volume(r,GNk,TGm1,TGm2,TGm3,m,nk,x_pole,thetaNumber)
%Pw=[]; 
Gw=[];Tm1=[];Tm2=[];Tm3=[];
%Pwn=[];
Gwn=[];TGmn1=[];TGmn2=[];TGmn3=[];    
for mk=1:length(x_pole)%累加所有模态-线变面*length(x_pole)
        %Pw{mk}=PNk{1,mk}*exp(sqrt(-1)*m(nk)*linspace(0,2*pi,360));%%(James-3.5.7)%nnk=nnk+1->验证可以执行parfor
        %PPw(:,:,mk)=Pw{1,mk}(r);
        Gw{mk}=GNk{1,mk}*exp(sqrt(-1)*m(nk)*linspace(0,2*pi,thetaNumber));%%(James-3.5.7)Reduced Green’s function at a particular frequency
        %GGw(:,:,mk)=Gw{1,mk}(r);
        Tm1{mk}=TGm1{1,mk}*exp(sqrt(-1)*m(nk)*(linspace(0,2*pi,thetaNumber)));%%(James-3.5.7)Reduced Green’s function at a particular frequency
        Tm2{mk}=TGm2{1,mk}*exp(sqrt(-1)*m(nk)*(linspace(0,2*pi,thetaNumber)));%%(James-3.5.7)Reduced Green’s function at a particular frequency
        Tm3{mk}=TGm3{1,mk}*exp(sqrt(-1)*m(nk)*(linspace(0,2*pi,thetaNumber)));%%(James-3.5.7)Reduced Green’s function at a particular frequency
        
        %TTm(:,:,mk)=Tm{1,mk}(r);
        %Pwn(mk)=max(sqrt(real(PNk{1,mk}(r)).^2));
%         Gwn(mk)=max(sqrt(real(GNk{1,mk}(r)).^2));
%         TGmn1(mk)=max(sqrt(real(TGm1{1,mk}(r)).^2));
%         TGmn2(mk)=max(sqrt(real(TGm2{1,mk}(r)).^2));
%         TGmn3(mk)=max(sqrt(real(TGm3{1,mk}(r)).^2));
        Gwn(mk)=max(abs(GNk{1,mk}(r)));
        TGmn1(mk)=max(abs(TGm1{1,mk}(r)));
        TGmn2(mk)=max(abs(TGm2{1,mk}(r)));
        TGmn3(mk)=max(abs(TGm3{1,mk}(r)));
    end
