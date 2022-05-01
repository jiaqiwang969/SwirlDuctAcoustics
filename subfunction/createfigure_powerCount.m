function createfigure_powerCount(XMatrix1, YMatrix1, ZMatrix1,color,Name,ztitle) %,
%CREATEFIGURE(XMatrix1, YMatrix1, ZMatrix1)
%  XMATRIX1:  x 数据的矩阵
%  YMATRIX1:  y 数据的矩阵
%  ZMATRIX1:  z 数据的矩阵

%  由 MATLAB 于 26-Oct-2018 11:22:07 自动生成

% 创建 figure
figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);

% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% 使用 plot3 的矩阵输入创建多行
plot31 = plot3(XMatrix1,YMatrix1,ZMatrix1,'Parent',axes1,'LineStyle','none',...
    'Color',[0 0 0]);
switch color
    case 1
set(plot31,'DisplayName',Name,'MarkerFaceColor',[1 0 0],...
    'MarkerSize',1,...
    'Marker','*',...
    'LineWidth',2);
    case 2
set(plot31,'DisplayName',Name,'MarkerFaceColor',[0 1 0],...
    'MarkerSize',2,...
    'Marker','+',...
    'LineWidth',2,...
    'Color',[1 0 0]);
    case 3
set(plot31,'DisplayName',Name,'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerSize',4,...
    'Marker','o');
end

% 创建 zlabel
zlabel({ztitle},'FontSize',11);

% 创建 ylabel
ylabel({'Imag(k)'},'FontSize',11);

% 创建 xlabel
xlabel({'Real(k)'},'FontSize',11);

% 取消以下行的注释以保留坐标区的 X 范围
xlim(axes1,[-40 40]);
% 取消以下行的注释以保留坐标区的 Y 范围
ylim(axes1,[-150 150]);
% 取消以下行的注释以保留坐标区的 Z 范围
zlim(axes1,[0 1]);
view(axes1,[-40.7000000000003 29.2]);
grid on;
