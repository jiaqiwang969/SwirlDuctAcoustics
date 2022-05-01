  % 保存图像
  function [h1,h2]=saveFun_RotorNoise(r,m,Tm,Gw,TGmn,Gwn,Tr,Omag,save_directory,Boundary,Type); 

    [rho,theta]=meshgrid(r',linspace(0,2*pi,360)); %生成极坐标网格
    [x,y]=pol2cart(theta,rho); %将极坐标网格转化为直角坐标网格

    
    h1=figure('Visible', 'on');
    subplot(2,2,1);contour(x,y,real(Gw)','Fill','on');axis equal
    colormap Parula;colorbar('location','EastOutside');
    title({['Plot of Real(Gw)'];[char(Type(Boundary)),'-Tr=',num2str(Tr),'-Omag=',num2str(Omag)]})
    subplot(2,2,2);contour(x,y,imag(Gw)','Fill','on');axis equal
    colormap Parula;colorbar('location','EastOutside');
    title({['Plot of Imag(Gw)'];[char(Type(Boundary)),'-Tr=',num2str(Tr),'-Omag=',num2str(Omag)]})
    subplot(2,2,[3,4]);
    bar(m,Gwn);set(gca,'YMinorTick','on','YScale','log');
    %set(gca,'XTick',m);
    set(gca,'Ygrid','on') 
    title({['Plot of maxGw for each azimuthal number n'];[char(Type(Boundary)),'-Tr=',num2str(Tr),'-Omag=',num2str(Omag)]},'FontSize',14)
    xlabel('Mode Number：m');ylabel('Amplitude');

    h2=figure('Visible', 'on');
    subplot(2,2,1);contour(x,y,real(Tm)','Fill','on');axis equal
    colormap Parula;colorbar('location','EastOutside');
    title({['Plot of Real(TGm)'];[char(Type(Boundary)),'-Tr=',num2str(Tr),'-Omag=',num2str(Omag)]})
    subplot(2,2,2);contour(x,y,imag(Tm)','Fill','on');axis equal
    colormap Parula;colorbar('location','EastOutside');
    title({['Plot of Imag(TGm)'];[char(Type(Boundary)),'-Tr=',num2str(Tr),'-Omag=',num2str(Omag)]})
    subplot(2,2,[3,4]);
    bar(m,TGmn);set(gca,'YMinorTick','on','YScale','log');
    %set(gca,'XTick',m);
    set(gca,'Ygrid','on') 
    title({['Plot of maxTGm for each azimuthal number n'];[char(Type(Boundary)),'-Tr=',num2str(Tr),'-Omag=',num2str(Omag)]},'FontSize',14)
    xlabel('Mode Number：m');ylabel('Amplitude');

    saveas(h1,[save_directory,'\','Gw',char(Type(Boundary)),'-Tr=',num2str(Tr),'-Omag=',num2str(Omag),'.png'])
    saveas(h2,[save_directory,'\','TGm',char(Type(Boundary)),'-Tr=',num2str(Tr),'-Omag=',num2str(Omag),'.png'])
    %saveas(h1,[save_directory,'\','Pw',char(Type(Boundary)),'-Tr=',num2str(Tr),'-Omag=',num2str(Omag),'.fig'])
    %saveas(h2,[save_directory,'\','Gw',char(Type(Boundary)),'-Tr=',num2str(Tr),'-Omag=',num2str(Omag),'.fig'])
  end