%%----------------------Pwer by SJTU. CopyLicense------------------------%%
%%fix soft boundary conditionbug,2018-10-7,wjq
%%add Logarithmic entropy condtion,2018-10-9,wjq
%%numerical green's function calculation,2018-10-10,wjq
%%colour plot of  Green¡¯s function Pw,2018-10-12,wjq
%%mode bar  at a particular frequency ,2018-10-13,wjq
%%Note: All dimensionless; applicable to Mx with r-varying shear flow
%%Derivation of GreenFunction
%%par-for is possible
%Generation of monopole or dipole stationary animation in the pipe
%%----------------------Pwer by SJTU. CopyLicense------------------------%%
clc;clear;close all
addpath(genpath(pwd));
warning('off')
% % //======Save=============== //
save_directory = ['result02-swirl/','output',date];
if ~exist(save_directory) mkdir(save_directory);else disp('file exit');end

r_pole=[0.8];%Pulse point position - upstream and downstream control
x_pole=linspace(-5,5,20);x_pole1=x_pole(find(x_pole<=0));x_pole2=x_pole(find(x_pole>0));
Entropy=0; %0:constant entropy condtion;1-2:logatithmic entropy condtion
Boundary=[1]; Type={'Hard Wall';'Lined Outer Wall';'Lined Inner Wall';'Lined Outer&Inner Wall'};
z_t=1-2*sqrt(-1);z_h=1-2*sqrt(-1);
beta=[0.3];
N =31;Ratio=0.6; [D,r] = cheb(N,Ratio,1);
l=4;probeNumber = 64; circshiftAngle = 6;thetaNumber=probeNumber*circshiftAngle;
w=25;Tr=0.2;Omag=0.4;Mx=0.4*ones(N+1,1);M_theta=[Tr./r(1:(N+1)/4*3);(Tr./r((N+1)/4*3)-Omag.*r((N+1)/4*3))+Omag.*r((N+1)/4*3+1:end)]; %
%Tr=0.0;Omag=0.5; Mx=0.001*ones(N+1,1);M_theta=Tr./r+Omag*r;
%When the spin coefficient is too large, then the effect must be calculated; Mx cannot be zero
[c02,rou0,P0,s0]=entropyPara(r,N,Ratio,Omag,Tr,Entropy,beta);
m=[-35:20];%52 53
tic
parfor nk=1:length(m)

    [V,lam]=eigfun_AB(r,D,N,w,m(nk),Ratio,Mx,M_theta,rou0,P0,c02,Boundary,z_t,z_h);%characteristics
    crLayer=[min((w-m(nk)*M_theta./r)./Mx);max((w-m(nk)*M_theta./r)./Mx)];
    if mean(abs(Mx))>0.02 cutOffLine=GMM_Cluster(lam,crLayer);  elseif  mean(abs(Mx))<0.02 cutOffLine=0;end
    [mode1]=eig_choose(V,N,lam,w,Mx,Tr,Omag,-2,40,-2,35,crLayer,cutOffLine,1); %for characteristics
    [mode2]=eig_choose(V,N,lam,w,Mx,Tr,Omag,-2,40,-2,35,crLayer,cutOffLine,2); %for characteristics
    [G_nm1,Tgm11,Tgm12,Tgm13]=greenfun_dipoleNoise(r,Boundary,m(nk),Ratio,w,Tr,Omag,Mx,c02,rou0,lam(mode1),z_t,z_h,r_pole,x_pole1,0,45,90);%green calculation£¨1*length(x_pole) cell£©
    [G_nm2,Tgm21,Tgm22,Tgm23]=greenfun_dipoleNoise(r,Boundary,m(nk),Ratio,w,Tr,Omag,Mx,c02,rou0,lam(mode2),z_t,z_h,r_pole,x_pole2,0,45,90);%green calculation£¨1*length(x_pole) cell£©
    [GNk1,TGm11,TGm12,TGm13]=cheb_cumKxCell(G_nm1,Tgm11,Tgm12,Tgm13,Ratio,length(x_pole1),length(mode1));%Sum all of wavenumbers*length(x_pole)
    [GNk2,TGm21,TGm22,TGm23]=cheb_cumKxCell(G_nm2,Tgm21,Tgm22,Tgm23,Ratio,length(x_pole2),length(mode2));%Sum all of wavenumbers*length(x_pole)
    [Gw1{1,nk},Tm11{1,nk},Tm12{1,nk},Tm13{1,nk},Gwn1{nk},TGmn11{nk},TGmn12{nk},TGmn13{nk}]=greenfun_span2volume(r,GNk1,TGm11,TGm12,TGm13,m,nk,x_pole1,thetaNumber);
    [Gw2{1,nk},Tm21{1,nk},Tm22{1,nk},Tm23{1,nk},Gwn2{nk},TGmn21{nk},TGmn22{nk},TGmn23{nk}]=greenfun_span2volume(r,GNk2,TGm21,TGm22,TGm23,m,nk,x_pole2,thetaNumber);

end
toc
%clearvars -except save_directory Boundary Type Tr Omag r m w Ratio x_pole1 x_pole2 x_pole ...
%Gw1 Tm11 Tm12 Tm13 Gwn1 TGmn11 TGmn12 TGmn13...
%Gw2 Tm21 Tm22 Tm23 Gwn2{nk} TGmn21 TGmn22 TGmn23...
%rou0 P0 s0 Mx M_theta



[GGw1,TTm11,TTm12,TTm13]=cheb_cumModeCell(r,Gw1,Tm11,Tm12,Tm13,Ratio,length(x_pole1));
[GGw2,TTm21,TTm22,TTm23]=cheb_cumModeCell(r,Gw2,Tm21,Tm22,Tm23,Ratio,length(x_pole2));
GGwn1=cell2mat(Gwn1');TTGmn11=cell2mat(TGmn11');TTGmn12=cell2mat(TGmn12');TTGmn13=cell2mat(TGmn13');
GGwn2=cell2mat(Gwn2');TTGmn21=cell2mat(TGmn21');TTGmn22=cell2mat(TGmn22');TTGmn23=cell2mat(TGmn23');
GGw=cat(3,GGw1,GGw2);TTm1=cat(3,TTm11,TTm21);TTm2=cat(3,TTm12,TTm22);TTm3=cat(3,TTm13,TTm23);
GGwn=cat(2,GGwn1,GGwn2);TTGmn1=cat(2,TTGmn11,TTGmn21);TTGmn2=cat(2,TTGmn12,TTGmn22);TTGmn3=cat(2,TTGmn13,TTGmn23);
[rou0_3D,P0_3D,s0_3D,Mx_3D,M_theta_3D]=meanflow_parameter(rou0,P0,s0,Mx,M_theta,x_pole,thetaNumber);


GGwArray=zeros(size(GGw));
TTm1Array=zeros(size(TTm1));
TTm2Array=zeros(size(TTm2));
TTm3Array=zeros(size(TTm3));

for k = 1:probeNumber
    GGwArray=GGwArray+circshift(GGw,circshiftAngle*(k-1),2)*exp(i*2*pi/probeNumber*k*l);
    TTm1Array=TTm1Array+circshift(TTm1,circshiftAngle*(k-1),2)*exp(i*2*pi/probeNumber*k*l);
    TTm2Array=TTm2Array+circshift(TTm2,circshiftAngle*(k-1),2)*exp(i*2*pi/probeNumber*k*l);
    TTm3Array=TTm3Array+circshift(TTm3,circshiftAngle*(k-1),2)*exp(i*2*pi/probeNumber*k*l);
end

tic

%[h1,h2]=saveFun_dipoleNoise(r,m,TTm(:,:,2),GGw(:,:,2),TTGmn(:,2),GGwn(:,2),Tr,Omag,save_directory,Boundary,Type);
%[h1,h2]=saveFun_dipoleNoise(r,m,TTm(:,:,end),GGw(:,:,end),TTGmn(:,end),GGwn(:,end),Tr,Omag,save_directory,Boundary,Type);

pltPlot_greenswirl(w,r,m,GGw,TTm1,TTm2,TTm3,GGwArray,TTm1Array,TTm2Array,TTm3Array,rou0_3D,P0_3D,s0_3D,Mx_3D,M_theta_3D,Tr,Omag,save_directory,Boundary,Type,x_pole,thetaNumber)
toc
%save([save_directory,'/',date,char(Type(Boundary)),'-Tr=',num2str(Tr),'-Omag=',num2str(Omag),'.mat'])
