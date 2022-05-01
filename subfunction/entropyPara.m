function [c02,rou0,P0,s0]=entropyPara(r,N,Ratio,Omag,Tr,Entropy,beta)     
Cp=1.005;Cv=0.718;gama=1.4;
switch  Entropy
    case {0} %Isentropic
        c02=1-1/2*(gama-1)*(Omag^2*(1-r.^2)-Tr^2*(1-1./(r.^2)))+2*(gama-1)*Omag*Tr*log(r);
        rou0=c02.^(1/(gama-1)); P0=1/gama*(c02.^(gama/(gama-1)));%Isentropic(非均熵多一个参数beita-3.36)
             
    case {1}  %by Guan 
        EntropyFun=chebop(Ratio,1);EntropyFun.rbc=1;EntropyFun.op=@(x,rou0) diff(rou0)-(beta/gama)./x*rou0-(Tr./x+Omag*x)^2*x^(beta-1)*rou0^(2-gama);
        rou0=EntropyFun\0;rou0=rou0(r);
        P0=1/gama*(rou0.^gama)./(r.^beta);   %no-Isentropic-3.36
        c02=gama*P0./rou0;
    case {2}  %by James   
        x = chebfun('x',[Ratio 1]);
        c02 = x^(-beta/Cp)-(gama-1)*x^(-beta/Cp)*cumsum(x^(beta/Cp)*(Tr/x+Omag*x)^2/x);c02=c02(r);%(James-2.1.19)
        rou0=(c02.*r.^(beta/Cv)).^(1/(gama-1));%(James-2.1.20)
        P0=(rou0.^gama)./gama.*exp(-log(r.^beta)./Cv); %(James-2.1.12) s0=-log(r^beta)
    case {3}  %for James  将gama->1 from Guan 参考熵为0=ln（1）
        EntropyFun=chebop(Ratio,1);EntropyFun.rbc=1;EntropyFun.op=@(x,rou0) diff(rou0)-(beta/gama)./x*rou0-1/gama*(Tr./x+Omag*x)^2*x^(beta-1)*rou0^(2-gama);%add 1/gama
        rou0=EntropyFun\0;rou0=rou0(r);
        P0=(rou0.^gama)./(r.^beta);   %no-Isentropic-3.36
        c02=gama*P0./rou0;
    case {4}  %for Guan2 same as Tam --3.87  rou0=1
        rou0=1*ones(N+1,1);
        P0=-(Omag^2/2-Tr^2/2-Omag^2*r.^2/2-2*Omag*Tr*log(r)+1/2*Tr^2./r.^2)+1/gama;  %(Guan-3.89)
        c02=-gama*(Omag^2/2-Tr^2/2-Omag^2*r.^2/2-2*Omag*Tr*log(r)+Tr^2/2/r.^2)+1;   %(Guan-3.90)
end
        s0=log(gama)+log(P0./(rou0.^gama));
end