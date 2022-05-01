function SpecialModeEigenValue(r,lam,N,V,n,r1,i1)
type={'脉动压力';'径向脉动速度';'周向脉动速度';'轴向脉动速度';'脉动压力';'径向涡量';'周向涡量';'轴向涡量';'熵值';};

figure('InvertHardcopy','off','Color',[1 1 1]);
specialMode= find(real(lam)>r1-2&real(lam)<r1+2&imag(lam)>i1-2&imag(lam)<i1+2);
[C,N]=min(lam(specialMode)-(r1+i*i1));
Mode=specialMode(N);
for kk=1:length(Mode)
plot(r,imag(V{1,Mode(kk)}(n,:)),'LineWidth',3);hold on;plot(r,real(V{1,Mode(kk)}(n,:)),'--','LineWidth',3);
end
title(num2str(lam(Mode)));
xlabel('r');
ylabel(type(n));