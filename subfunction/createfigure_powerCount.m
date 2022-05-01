function createfigure_powerCount(XMatrix1, YMatrix1, ZMatrix1,color,Name,ztitle) %,
%CREATEFIGURE(XMatrix1, YMatrix1, ZMatrix1)
%  XMATRIX1:  x ���ݵľ���
%  YMATRIX1:  y ���ݵľ���
%  ZMATRIX1:  z ���ݵľ���

%  �� MATLAB �� 26-Oct-2018 11:22:07 �Զ�����

% ���� figure
figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);

% ���� axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% ʹ�� plot3 �ľ������봴������
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

% ���� zlabel
zlabel({ztitle},'FontSize',11);

% ���� ylabel
ylabel({'Imag(k)'},'FontSize',11);

% ���� xlabel
xlabel({'Real(k)'},'FontSize',11);

% ȡ�������е�ע���Ա����������� X ��Χ
xlim(axes1,[-40 40]);
% ȡ�������е�ע���Ա����������� Y ��Χ
ylim(axes1,[-150 150]);
% ȡ�������е�ע���Ա����������� Z ��Χ
zlim(axes1,[0 1]);
view(axes1,[-40.7000000000003 29.2]);
grid on;
