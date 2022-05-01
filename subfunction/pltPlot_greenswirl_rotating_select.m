function pltPlot_greenswirl_rotating(r,t,v1,save_directory,Boundary,Type,x_pole,name)
tic
for solutime=1:length(t)
%% 生成网格矩阵
    [theta,rho,z]=meshgrid(linspace(0,2*pi,180),r',x_pole); %生成极坐标网格
    [y,x,z]=pol2cart(theta,rho,z); %将极坐标网格转化为直角坐标网格  
% 导入数据，对网格赋值
 %%
     tsignal.cubes(solutime).zonename='mysurface zone';
     tsignal.cubes(solutime).x=x;    %size 3x3 
     tsignal.cubes(solutime).y=y;    %size 3x3
     tsignal.cubes(solutime).z=z;    %size 3x3
      tsignal.cubes(solutime).v(1,:,:,:)=v1{solutime};%可以为三维输入，p只有一维
     tsignal.cubes(solutime).solutiontime=solutime;
end
 %% 整合数据，生成文件
%title=''; 

NAME = [strrep(strrep(char(datetime),':','-'),' ','-'),'Ro-GreenSwirl',name];  %存储文件夹   
output_file_name=[save_directory,'\',NAME,'.plt']; 
tsignal.Nvar=4;     
tsignal.varnames={'x','y','z',NAME};
mat2tecplot(tsignal,output_file_name);
toc
end
 
 
