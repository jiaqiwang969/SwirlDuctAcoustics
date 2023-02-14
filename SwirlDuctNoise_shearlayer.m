clc;
clear;
% close all

addpath(genpath(pwd));
warning('off')
% % //======Save=============== //
save_directory = ['result02-swirl/','output',date];
if ~exist(save_directory) mkdir(save_directory);else disp('file exit');end

%%
r_pole=[0.8];%Pulse point position - upstream and downstream control
x_pole=linspace(-2,2,100);x_pole1=x_pole(find(x_pole<=0));x_pole2=x_pole(find(x_pole>0));
Entropy=0; %0:constant entropy condtion;1-2:logatithmic entropy condtion
Boundary=[1]; Type={'Hard Wall';'Lined Outer Wall';'Lined Inner Wall';'Lined Outer&Inner Wall'};
z_t=1-2*sqrt(-1);z_h=1-2*sqrt(-1);
beta=[0.3];
N =31;Ratio=0.6; [D,r] = cheb(N,Ratio,1);
l=4;probeNumber = 64; circshiftAngle = 6;thetaNumber=probeNumber*circshiftAngle;
w=25;Tr=0.1;Omag=0.1;M_theta=Tr./r+Omag*r; %
%%%%
%Mx通过插值得到
%Mx=0.4*ones(N+1,1);
load fit_3.mat
xr0=xr;u00=u0;
[fitresult0, gof] = createFit(xr, u0)
Mx0=fitresult0(r);
load fit_1.mat
xr1=xr;u01=u0;
[fitresult1, gof1] = createFit(xr, u0)
Mx1=fitresult1(r);

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xr, u0 );
% legend( h, 'u0 vs. r', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'r', 'Interpreter', 'none' );
% ylabel( 'u0', 'Interpreter', 'none' );
% grid on

% syms t
% % f=fittype('a*exp(-b*t)*sin((2*pi/c)*t*exp(-d*t)+e)+f','independent','t','coefficients',{'a','b','c','d','e','f'});
% f=fittype('a./(b+c*exp(d*t-e))+f','independent','t','coefficients',{'a','b','c','d','e','f'});
% cfun=fit(xr,u0,f,'Normalize','on','Robust', 'Bisquare');
% 


[c02,rou0,P0,s0]=entropyPara(r,N,Ratio,Omag,Tr,Entropy,beta);

m=[-35:35];%52 53
Nk=50

for nk=1:length(m)
%     figure
%     c=colormap(jet(Nk+1));
    for k=1:Nk+1
        u0=u00+(u01-u00)/Nk*(k-1);
        [fitresult, gof] = createFit(xr, u0);
        Mx=smooth(fitresult(r),0.1);
        [V(k,nk,:),lam(:,nk,k)]=eigfun_AB(r,D,N,w,m(nk),Ratio,Mx,M_theta,rou0,P0,c02,Boundary,z_t,z_h);%characteristics
%         subplot(1,2,1)
%         plot(r,Mx,'Color',c(k,:));hold on
%         xlabel("r");ylabel("u0");
%         subplot(1,2,2)
%         plot(real(lam(:,nk,k)),imag(lam(:,nk,k)),'.','Color',c(k,:));hold on
%         xlim([-150,60])
%         ylim([-50,50])
%         xlabel("real(k)");ylabel("imag(k)");
    end
end

dr=0.5;


for nk=1:length(m)
    aline_up=find(real(lam(:,nk,1))+2*imag(lam(:,nk,1))+10<0 & abs(imag(lam(:,nk,1)))<50  & abs(real(lam(:,nk,1)))<1000);
    aline_down=find(real(lam(:,nk,1))+2*imag(lam(:,nk,1))+10>0 & abs(imag(lam(:,nk,1)))<50 & abs(real(lam(:,nk,1)))<1000);

%     figure
%     c=colormap(jet(Nk+1));
    for k=1:Nk
    %在每个特征值lam(:,mk,k)上画一个半价为r的小圆圈，若后面lam(:,mk,k+1)包含在里面，且个数为1，则保留！
    relation_up=abs(repmat(lam(aline_up,nk,k).',length(lam(:,nk,k+1)),1) - repmat(lam(:,nk,k+1),1,length(aline_up)))<dr;
    relation_down=abs(repmat(lam(aline_down,nk,k).',length(lam(:,nk,k+1)),1) - repmat(lam(:,nk,k+1),1,length(aline_down)))<dr;
    relation_sum_up=sum(relation_up,2); %在lam(:,mk,k+1)的小圈圈内，包含lam(:,mk,k)的个数，1个说明lam(:,mk,k+1)对应的是acoustic mode
    relation_sum_down=sum(relation_down,2); %在lam(:,mk,k+1)的小圈圈内，包含lam(:,mk,k)的个数，1个说明lam(:,mk,k+1)对应的是acoustic mode
    [aline_up]=find(relation_sum_up==1);
    [aline_down]=find(relation_sum_down==1);

%         subplot(1,2,2)
%         plot(real(lam(aline_up,nk,k+1)),imag(lam(aline_up,nk,k+1)),'.','Color',c(k,:));hold on
%         plot(real(lam(aline_down,nk,k+1)),imag(lam(aline_down,nk,k+1)),'x','Color',c(k,:));hold on
%         xlim([-150,60])
%         ylim([-50,50])
%         xlabel("real(k)");ylabel("imag(k)");
    end
    lam(aline_up,nk,k+1)
    lam(aline_down,nk,k+1)

    nk
    [G_nm1,Tgm11,Tgm12,Tgm13]=greenfun_dipoleNoise(r,Boundary,m(nk),Ratio,w,Tr,Omag,Mx,c02,rou0,lam(aline_up,nk,k+1),z_t,z_h,r_pole,x_pole1,0,45,90);%green calculation（1*length(x_pole) cell）
    [G_nm2,Tgm21,Tgm22,Tgm23]=greenfun_dipoleNoise(r,Boundary,m(nk),Ratio,w,Tr,Omag,Mx,c02,rou0,lam(aline_down,nk,k+1)  ,z_t,z_h,r_pole,x_pole2,0,45,90);%green calculation（1*length(x_pole) cell）
    [GNk1,TGm11,TGm12,TGm13]=cheb_cumKxCell(G_nm1,Tgm11,Tgm12,Tgm13,Ratio,length(x_pole1),length(lam(aline_up,nk,k+1)));%Sum all of wavenumbers*length(x_pole)
    [GNk2,TGm21,TGm22,TGm23]=cheb_cumKxCell(G_nm2,Tgm21,Tgm22,Tgm23,Ratio,length(x_pole2),length(lam(aline_down,nk,k+1)));%Sum all of wavenumbers*length(x_pole)
    [Gw1{1,nk},Tm11{1,nk},Tm12{1,nk},Tm13{1,nk},Gwn1{nk},TGmn11{nk},TGmn12{nk},TGmn13{nk}]=greenfun_span2volume(r,GNk1,TGm11,TGm12,TGm13,m,nk,x_pole1,thetaNumber);
    [Gw2{1,nk},Tm21{1,nk},Tm22{1,nk},Tm23{1,nk},Gwn2{nk},TGmn21{nk},TGmn22{nk},TGmn23{nk}]=greenfun_span2volume(r,GNk2,TGm21,TGm22,TGm23,m,nk,x_pole2,thetaNumber);

end




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


figure;

subplot(1,3,1)
bar(m,TTGmn1(:,1))
xlabel('Mode number')
grid minor
subplot(1,3,2)
bar(m,TTGmn2(:,1))
xlabel('Mode number')

grid minor
subplot(1,3,3)
bar(m,TTGmn3(:,1))
xlabel('Mode number')

grid minor






