function [mode]=eig_choose(V,N,lam,w,Mx,Tr,Omag,re0,re1,im0,im1,crLayer,cutOffLine,isdownStream);  %10 36 -2 35
%选特征:如果有旋流，则存在critical layer
if isdownStream==1
    cmode=find(real(lam)<crLayer(2)+5&real(lam)>crLayer(1)-5&imag(lam)<2&imag(lam)>-2);
    cutOnmode= find(real(lam)>(cutOffLine+re0)&real(lam)<(cutOffLine+re1)&imag(lam)>-0.5&imag(lam)<0.5&(real(lam).^2+imag(lam).^2)>1E-3);
    cutOffmode= find(real(lam)>(cutOffLine-2)&real(lam)<(cutOffLine+2)&imag(lam)>im0&imag(lam)<im1&(real(lam).^2+imag(lam).^2)>1E-3);
    mode=[cutOnmode;cutOffmode];
elseif isdownStream==-1
    cmode=find(real(lam)<crLayer(2)+5&real(lam)>crLayer(1)-5&imag(lam)<2&imag(lam)>-2);
    cutOnmode= find(real(lam)>(cutOffLine-re1)&real(lam)<(cutOffLine-re0)&imag(lam)>-0.5&imag(lam)<0.5&(real(lam).^2+imag(lam).^2)>1E-3);
    cutOffmode= find(real(lam)>(cutOffLine-2)&real(lam)<(cutOffLine+2)&imag(lam)>-im1&imag(lam)<-im0&(real(lam).^2+imag(lam).^2)>1E-3);
    mode=[cutOnmode;cutOffmode];
end

    for ii=1:length(cmode)
    mode(mode==cmode(ii))=[];
    end   
    Kx=lam(mode)
 end

% s = chebfun('s',[0 1]);
% z = join(c(1)+s*(c(2)-c(1)),c(2)+s*(c(3)-c(2)) , ...
% c(3)+s*(c(4)-c(3)), c(4)+s*(c(1)-c(4)));

% for kk=1:5*(N+1)
% vf{kk} = [reshape(V(:,kk),[5,N+1])];
% mvalue(kk,:)=meanvalue(vf');
% Feature=features_extract(vf);
% for ii=1:length(Feature)
%     feature(:,ii)=Feature{ii}(:,1);
% end
% figure
% plot(1:510,feature)
% end
% mvalue=meanvalue(vf')
% 
%     figure;plot(real(lam),imag(lam),'.');hold on
%    % title(num2str(m(nk)));
%     grid on;  
%     set(gca, 'XLim',[-100 200]); 
%     set(gca, 'YLim',[-100 100]);
%     xlabel Real, ylabel Imaginary    
