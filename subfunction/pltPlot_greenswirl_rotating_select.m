function pltPlot_greenswirl_rotating(r,t,v1,save_directory,Boundary,Type,x_pole,name)
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
     tsignal.cubes(solutime).solutiontime=solutime;
end
 %% �������ݣ������ļ�
%title=''; 

NAME = [strrep(strrep(char(datetime),':','-'),' ','-'),'Ro-GreenSwirl',name];  %�洢�ļ���   
output_file_name=[save_directory,'\',NAME,'.plt']; 
tsignal.Nvar=4;     
tsignal.varnames={'x','y','z',NAME};
mat2tecplot(tsignal,output_file_name);
toc
end
 
 
