function pltPlot_greenswirl(w,r,m,GGw,TTm1,TTm2,TTm3,GGwArray,TTm1Array,TTm2Array,TTm3Array,rou0_3D,P0_3D,s0_3D,Mx_3D,M_theta_3D,Tr,Omag,save_directory,Boundary,Type,x_pole,thetaNumber)
T1=8*pi/w;
t=linspace(0,T1,50);

for solutime=1:length(t)
%% �����������
    [theta,rho,z]=meshgrid(linspace(0,2*pi,thetaNumber),r',x_pole); %���ɼ���������
    [y,x,z]=pol2cart(theta,rho,z); %������������ת��Ϊֱ����������
    
% �������ݣ�������ֵ

 %%
     tsignal.cubes(solutime).zonename='mysurface zone';
     tsignal.cubes(solutime).x=x;    %size 3x3 
     tsignal.cubes(solutime).y=y;    %size 3x3
     tsignal.cubes(solutime).z=z;    %size 3x3
     tsignal.cubes(solutime).v(1,:,:,:)=real(GGw*exp(-sqrt(-1)*w*t(solutime)));%����Ϊ��ά���룬pֻ��һά
     tsignal.cubes(solutime).v(2,:,:,:)=real(TTm1*exp(-sqrt(-1)*w*t(solutime)));%0��
     tsignal.cubes(solutime).v(3,:,:,:)=real(TTm2*exp(-sqrt(-1)*w*t(solutime)));%45��
     tsignal.cubes(solutime).v(4,:,:,:)=real(TTm3*exp(-sqrt(-1)*w*t(solutime)));%90��
     tsignal.cubes(solutime).v(5,:,:,:)=rou0_3D;%ƽ�����ܶ�
     tsignal.cubes(solutime).v(6,:,:,:)=P0_3D;%ƽ����ѹ��
     tsignal.cubes(solutime).v(7,:,:,:)=s0_3D;%ƽ������ֵ
     tsignal.cubes(solutime).v(8,:,:,:)=Mx_3D;%ƽ���������ٶ�
     tsignal.cubes(solutime).v(9,:,:,:)=M_theta_3D;%ƽ�������ٶ�
     tsignal.cubes(solutime).v(10,:,:,:)=real(GGwArray*exp(-sqrt(-1)*w*t(solutime)));%����Ϊ��ά���룬pֻ��һά
     tsignal.cubes(solutime).v(11,:,:,:)=real(TTm1Array*exp(-sqrt(-1)*w*t(solutime)));%0��
     tsignal.cubes(solutime).v(12,:,:,:)=real(TTm2Array*exp(-sqrt(-1)*w*t(solutime)));%45��
     tsignal.cubes(solutime).v(13,:,:,:)=real(TTm3Array*exp(-sqrt(-1)*w*t(solutime)));%90��
     %tsignal.cubes(solutime).v(10,:,:,:)=Mr_3D;%ƽ���������ٶ�
     tsignal.cubes(solutime).solutiontime=solutime;
end
 %% �������ݣ������ļ�
%title=''; 
NAME = [date,'GreenSwirl',char(Type(Boundary)),'-Tr=',num2str(Tr),'-Omag=',num2str(Omag),'-',num2str(x_pole(1)),'-',num2str(x_pole(end))];  %�洢�ļ���   
NAME1 = ['GreenSwirl-Gw-0','-',date];  %�洢�ļ���   
NAME2 = ['GreenSwirl-Tm-0','-',date];  %�洢�ļ���   
NAME3 = ['GreenSwirl-Tm-45','-',date];  %�洢�ļ���   
NAME4 = ['GreenSwirl-Tm-90','-',date];  %�洢�ļ���  
NAME5 = ['GreenSwirl-rou0','-',date];  %�洢�ļ���   
NAME6 = ['GreenSwirl-P0_3D','-',date];  %�洢�ļ���   
NAME7 = ['GreenSwirl-s0_3D','-',date];  %�洢�ļ���   
NAME8 = ['GreenSwirl-Mx_3D','-',date];  %�洢�ļ���   
NAME9 = ['GreenSwirl-M_theta_3D','-',date];  %�洢�ļ���   
NAME10 =['GreenSwirl-GwmArray','-',date];
NAME11 =['GreenSwirl-Tm1Array','-',date];
NAME12 =['GreenSwirl-Tm2Array','-',date];
NAME13 =['GreenSwirl-Tm3Array','-',date];
%NAME10 = ['GreenSwirl-Mr_3D','-',date];  %�洢�ļ���   �0�8


output_file_name=[save_directory,'/',NAME,'.plt']; 
tsignal.Nvar=16;     
tsignal.varnames={'x','y','z',NAME1,NAME2,NAME3,NAME4,NAME5,NAME6,NAME7,NAME8,NAME9,NAME10,NAME11,NAME12,NAME13};
mat2tecplot(tsignal,output_file_name);

end
 
 
