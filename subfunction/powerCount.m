function [Hp,Hw,Hs,V]=powerCount(r,N,D,V,lam,m,rou0,P0,color,displayName)
%code from(Guan-3.96-3.101)
Cp=1.005;Cv=0.718;
for kk=1:length(lam)
V{kk}(6,:)=V{kk}(3,:).'./r+D*V{kk}(3,:).'-V{kk}(2,:).'./r*i*m;
V{kk}(7,:)=V{kk}(4,:).'./r*i*m-V{kk}(3,:).'.*lam(kk)*i;  
V{kk}(8,:)=V{kk}(2,:).'.*i*lam(kk)-D*V{kk}(4,:).';
V{kk}(9,:)=Cv*V{kk}(1,:).'./P0-Cp*V{kk}(5,:).'./rou0;
v1=abs(V{kk}(1,:)).';%脉动压强 
v2=abs(V{kk}(6,:)).';%x方向涡量
v3=abs(V{kk}(7,:)).';%%r方向涡量
v4=abs(V{kk}(8,:)).';%%theta方向涡量
v5=abs(V{kk}(9,:)).';%%脉动熵
V1=sum(v1);
V2=sum(v2+v3+v4);
V3=sum(v5);
Hp(kk,1)=V1/(V1+V2+V3);
Hw(kk,1)=V2/(V1+V2+V3);
Hs(kk,1)=V3/(V1+V2+V3);
end
createfigure_powerCount(real(lam),imag(lam),Hp,color,displayName,'P_P_r_e_s_s_u_r_e');
createfigure_powerCount(real(lam),imag(lam),Hw,color,displayName,'P_V_o_r_t_i_c_i_t_y');
createfigure_powerCount(real(lam),imag(lam),Hs,color,displayName,'P_E_n_t_r_o_p_y');
createfigure_powerCount(real(lam),imag(lam),Hs,color,displayName,'Power_P_r_e_s_s_u_r_e');

    
end
