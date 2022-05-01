function [mode,len]=eig_choose(V,N,lam,w,Mx,Tr,Omag,re0,re1,im0,im1,crLayer,cutOffLine,n_max,isdownStream);  %10 36 -2 35
%选特征:如果有旋流，则存在critical layer
len=[];
if isdownStream==2
    cmode=find(real(lam)<crLayer(2)+5&real(lam)>crLayer(1)-5&imag(lam)<2&imag(lam)>-2);
    cutOnmode= find(real(lam)>(cutOffLine+re0)&real(lam)<(cutOffLine+re1)&imag(lam)>-0.5&imag(lam)<0.5&(real(lam).^2+imag(lam).^2)>1E-3);
    cutOffmode= find(real(lam)>(cutOffLine-2)&real(lam)<(cutOffLine+2)&imag(lam)>im0&imag(lam)<im1&(real(lam).^2+imag(lam).^2)>1E-3);
    mode=[cutOnmode;cutOffmode];

elseif isdownStream==1
    cmode=find(real(lam)<crLayer(2)+5&real(lam)>crLayer(1)-5&imag(lam)<2&imag(lam)>-2);
    cutOnmode= find(real(lam)>(cutOffLine-re1)&real(lam)<(cutOffLine-re0)&imag(lam)>-0.5&imag(lam)<0.5&(real(lam).^2+imag(lam).^2)>1E-3);
    cutOffmode= find(real(lam)>(cutOffLine-2)&real(lam)<(cutOffLine+2)&imag(lam)>-im1&imag(lam)<-im0&(real(lam).^2+imag(lam).^2)>1E-3);
    mode=[cutOnmode;cutOffmode];
    
end

    for ii=1:length(cmode)
    mode(mode==cmode(ii))=[];
    end
    for kk=1:length(mode)
       len(kk)=length(roots(chebfun(imag(V{1,mode(kk)}(1,:)'))));
    end   
      mode(find((len>n_max)==1))=[];

    Kx=lam(mode)
 end

