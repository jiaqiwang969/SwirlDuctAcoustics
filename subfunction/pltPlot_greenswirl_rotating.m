function pltPlot_greenswirl_rotating(t,r,m,v1,v2,v3,v4,rou0_3D,P0_3D,s0_3D,Mx_3D,M_theta_3D,Tr,Omag,save_directory,Boundary,Type,x_pole)
tic
for solutime=1:length(t)
%% �����������
    [theta,rho,z]=meshgrid(linspace(0,2*pi,180),r',x_pole); %���ɼ���������
    [y,x,z]=pol2cart(theta,rho,z); %������������ת��Ϊֱ����������
    
% �������ݣ�������ֵ

 %%
     tsignal.cubes(solutime).zonename='mysurface zone';
     tsignal.cubes(solutime).x=x;    %size 3x3 
     tsignal.cubes(solutime).y=y;    %size 3x3
     tsignal.cubes(solutime).z=z;    %size 3x3
      tsignal.cubes(solutime).v(1,:,:,:)=v1{solutime};%����Ϊ��ά���룬pֻ��һά
      tsignal.cubes(solutime).v(2,:,:,:)=v2{solutime};%0��
      tsignal.cubes(solutime).v(3,:,:,:)=v3{solutime};%45��
      tsignal.cubes(solutime).v(4,:,:,:)=v4{solutime};%90��
     tsignal.cubes(solutime).v(5,:,:,:)=rou0_3D;%ƽ�����ܶ�
     tsignal.cubes(solutime).v(6,:,:,:)=P0_3D;%ƽ����ѹ��
     tsignal.cubes(solutime).v(7,:,:,:)=s0_3D;%ƽ������ֵ
     tsignal.cubes(solutime).v(8,:,:,:)=Mx_3D;%ƽ���������ٶ�
     tsignal.cubes(solutime).v(9,:,:,:)=M_theta_3D;%ƽ�������ٶ�
     %tsignal.cubes(solutime).v(10,:,:,:)=Mr_3D;%ƽ���������ٶ�
     tsignal.cubes(solutime).solutiontime=solutime;
end
 %% �������ݣ������ļ�
%title=''; 
NAME = [date,'2Ro-GreenSwirl',char(Type(Boundary)),'-Tr=',num2str(Tr),'-Omag=',num2str(Omag),'-',num2str(x_pole(1)),'-',num2str(x_pole(end))];  %�洢�ļ���   
NAME1 = ['Ro-GreenSwirl-Gw-0','-',date];  %�洢�ļ���   
NAME2 = ['Ro-GreenSwirl-Tm-0','-',date];  %�洢�ļ���   
NAME3 = ['Ro-GreenSwirl-Tm-45','-',date];  %�洢�ļ���   
NAME4 = ['Ro-GreenSwirl-Tm-90','-',date];  %�洢�ļ���   
NAME5 = ['Ro-GreenSwirl-rou0','-',date];  %�洢�ļ���   
NAME6 = ['Ro-GreenSwirl-P0_3D','-',date];  %�洢�ļ���   
NAME7 = ['Ro-GreenSwirl-s0_3D','-',date];  %�洢�ļ���   
NAME8 = ['Ro-GreenSwirl-Mx_3D','-',date];  %�洢�ļ���   
NAME9 = ['Ro-GreenSwirl-M_theta_3D','-',date];  %�洢�ļ���   
%NAME10 = ['GreenSwirl-Mr_3D','-',date];  %�洢�ļ���   


output_file_name=[save_directory,'\',NAME,'.plt']; 
tsignal.Nvar=12;     
tsignal.varnames={'x','y','z',NAME1,NAME2,NAME3,NAME4,NAME5,NAME6,NAME7,NAME8,NAME9};
mat2tecplot(tsignal,output_file_name);
toc
end
 
 
