% Author: Weicen Dong <orange233@sjtu.edu.cn>
% Created Date: 2025/1/21

load('data_phase_QMT.mat')
col=C_vec;
col(:,(end-1)/2+1)=zeros(nc,1);col(:,(end-1)/2)=zeros(nc,1);col(:,(end-1)/2+2)=zeros(nc,1);
load('b_w_r.mat')
xte=-0.15; yte=1;
figure('Color','white'),h1=subplot(1,5,1);
pcolor(c_vec,d_vec,real(col'));xlabel('$n-m$',  'Interpreter', 'latex');ylabel('$b-a$',  'Interpreter', 'latex');
patch([-1 NaN],[0 NaN],[-2 NaN]),patch([-1 NaN],[0 NaN],[2 NaN]),xlim([-0.1 0.1])
shading flat;
text(0,0.05 ,'$a<b, C^{1}=-1$',  'Interpreter', 'latex','HorizontalAlignment', 'center', 'FontSize', 14);
text(0,-0.05 ,'$a>b, C^{1}=1$',  'Interpreter', 'latex','HorizontalAlignment', 'center', 'FontSize', 14);
xticks(-0.1:0.1:0.1);yticks(-0.1:0.1:0.1);box on
%colormap(h1,parula)
load('gg_sun.mat'),colormap(h1,gg_sun)
text(xte, yte, '(a)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');
xlim([-0.101,0.1]),box on

figure('Color','white')
t = tiledlayout(2, 5, 'TileSpacing', 'compact', 'Padding', 'loose');%TileSpacing: 控制子图之间的间距，可以是 'none', 'compact', 或 'loose'。•	Padding: 控制整个子图组的外边距。
% 子图编号:
% 1   2   3   4
% 5   6   7   8

% t.Children 顺序:
% 8   7   6   5
% 4   3   2   1
for i = 1:10
    nexttile;
    plot(rand(10,1)); % 示例内容
    title(['Subplot ' num2str(i)]);
end
axes_list = flipud(t.Children);% 翻转顺序，使其与创建顺序一致
xte=0.2; yte=0.9;

t1=1;a=0.1;c=0.2;b=0.2;d=0.1;
[C,k11,w,B] = Chern_mol(t1,a,b,c,d,L);C
[K1,K2]=meshgrid(k11,k11);
KX=K1;KY=(K1+2*K2)/sqrt(3);
h2 = axes_list(1); 
axes(h2); % 将其设为当前激活的子图
pcolor(KX(1:end-1,1:end-1),KY(1:end-1,1:end-1),real(B)),hold on
patch([-1.1 NaN],[1.1 NaN],[-4 NaN],'EdgeColor','interp'),patch([-1.1 NaN],[1.1 NaN],[4 NaN],'EdgeColor','interp')
shading interp,xlabel('k_1'),ylabel('k_2'),axis equal,axis off,colormap(h2,b_w_r)
text(0.5,0.9 ,'Re(B)','Units', 'normalized','HorizontalAlignment', 'center', 'FontSize', 10);
plot(KX(1,:),KY(1,:),'k'),plot(KX(end,:),KY(end,:),'k'),plot(KX(:,1),KY(:,1),'k'),plot(KX(:,end-1),KY(:,end-1),'k')
text(0.6,1.1 ,'$m<n, a>b$','Units', 'normalized','HorizontalAlignment', 'center', 'FontSize', 14,'Color','r',  'Interpreter', 'latex');
text(0.55,0.1 ,'k_1','Units', 'normalized', 'FontSize', 10);text(0.1,0.35 ,'k_2','Units', 'normalized', 'FontSize', 10);
text(KX(1,1)-0.3,KY(1,1)-0.3,'0', 'FontSize', 10,'HorizontalAlignment', 'center')
text(KX(1,end-1)-0.1,KY(1,end-1)-0.5,'1', 'FontSize', 10,'HorizontalAlignment', 'center')
text(KX(end-1,1)-0.3,KY(end-1,1)-0.3,'1', 'FontSize', 10,'HorizontalAlignment', 'center')
text(xte, 1.1, '(b)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');

h3 = axes_list(6);axes(h3);   % 将其设为当前激活的子图
pcolor(KX(1:end-1,1:end-1),KY(1:end-1,1:end-1),imag(B)),hold on
patch([-1.1 NaN],[1.1 NaN],[-4 NaN],'EdgeColor','interp'),patch([-1.1 NaN],[1.1 NaN],[4 NaN],'EdgeColor','interp')
shading interp,xlabel('k_1'),ylabel('k_2'),axis equal,axis off,colormap(h3,b_w_r)
plot(KX(1,:),KY(1,:),'k'),plot(KX(end,:),KY(end,:),'k'),plot(KX(:,1),KY(:,1),'k'),plot(KX(:,end-1),KY(:,end-1),'k')
text(0.5,0.9 ,'Im(B)','Units', 'normalized','HorizontalAlignment', 'center', 'FontSize', 10);
text(0.55,0.1 ,'k_1','Units', 'normalized', 'FontSize', 10);text(0.1,0.35 ,'k_2','Units', 'normalized', 'FontSize', 10);
text(KX(1,1)-0.3,KY(1,1)-0.3,'0', 'FontSize', 10,'HorizontalAlignment', 'center')
text(KX(1,end-1)-0.1,KY(1,end-1)-0.5,'1', 'FontSize', 10,'HorizontalAlignment', 'center')
text(KX(end-1,1)-0.3,KY(end-1,1)-0.3,'1', 'FontSize', 10,'HorizontalAlignment', 'center')
%text(xte, yte, '(f)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');

t1=1;a=0.2;c=0.1;b=0.2;d=0.1;
[C,k11,w,B] = Chern_mol(t1,a,b,c,d,L);C
h4=axes_list(2);axes(h4);   % 将其设为当前激活的子图
pcolor(KX(1:end-1,1:end-1),KY(1:end-1,1:end-1),real(B)),hold on
patch([-1.1 NaN],[1.1 NaN],[-4 NaN],'EdgeColor','interp'),patch([-1.1 NaN],[1.1 NaN],[4 NaN],'EdgeColor','interp')
shading interp,xlabel('k_1'),ylabel('k_2'),axis equal,axis off,colormap(h4,b_w_r)
plot(KX(1,:),KY(1,:),'k'),plot(KX(end,:),KY(end,:),'k'),plot(KX(:,1),KY(:,1),'k'),plot(KX(:,end-1),KY(:,end-1),'k')
text(0.5,0.9 ,'Re(B)','Units', 'normalized','HorizontalAlignment', 'center', 'FontSize', 10);
text(0.6,1.1 ,'$m>n, a>b$','Units', 'normalized','HorizontalAlignment', 'center', 'FontSize', 14,'Color','r',  'Interpreter', 'latex');
text(0.55,0.1 ,'k_1','Units', 'normalized', 'FontSize', 10);text(0.1,0.35 ,'k_2','Units', 'normalized', 'FontSize', 10);
text(KX(1,1)-0.3,KY(1,1)-0.3,'0', 'FontSize', 10,'HorizontalAlignment', 'center')
text(KX(1,end-1)-0.1,KY(1,end-1)-0.4,'1', 'FontSize', 10,'HorizontalAlignment', 'center')
text(KX(end-1,1)-0.3,KY(end-1,1)-0.3,'1', 'FontSize', 10,'HorizontalAlignment', 'center')
text(xte, 1.1, '(c)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');

h5=axes_list(7);axes(h5);   % 将其设为当前激活的子图
pcolor(KX(1:end-1,1:end-1),KY(1:end-1,1:end-1),imag(B)),hold on
patch([-1.1 NaN],[1.1 NaN],[-4 NaN],'EdgeColor','interp'),patch([-1.1 NaN],[1.1 NaN],[4 NaN],'EdgeColor','interp')
shading interp,xlabel('k_1'),ylabel('k_2'),axis equal,axis off,colormap(h5,b_w_r)
plot(KX(1,:),KY(1,:),'k'),plot(KX(end,:),KY(end,:),'k'),plot(KX(:,1),KY(:,1),'k'),plot(KX(:,end-1),KY(:,end-1),'k')
text(0.5,0.9 ,'Im(B)','Units', 'normalized','HorizontalAlignment', 'center', 'FontSize', 10);
text(0.55,0.1 ,'k_1','Units', 'normalized', 'FontSize', 10);text(0.1,0.35 ,'k_2','Units', 'normalized', 'FontSize', 10);
text(KX(1,1)-0.3,KY(1,1)-0.3,'0', 'FontSize', 10,'HorizontalAlignment', 'center')
text(KX(1,end-1)-0.1,KY(1,end-1)-0.4,'1', 'FontSize', 10,'HorizontalAlignment', 'center')
text(KX(end-1,1)-0.3,KY(end-1,1)-0.3,'1', 'FontSize', 10,'HorizontalAlignment', 'center')

t1=1;a=0.1;c=0.2;b=0.1;d=0.2;
[C,k11,w,B] = Chern_mol(t1,a,b,c,d,L);C
h6=axes_list(3);axes(h6);   % 将其设为当前激活的子图
pcolor(KX(1:end-1,1:end-1),KY(1:end-1,1:end-1),real(B)),hold on
patch([-1.1 NaN],[1.1 NaN],[-4 NaN],'EdgeColor','interp'),patch([-1.1 NaN],[1.1 NaN],[4 NaN],'EdgeColor','interp')
shading interp,xlabel('k_1'),ylabel('k_2'),axis equal,axis off,colormap(h6,b_w_r)
plot(KX(1,:),KY(1,:),'k'),plot(KX(end,:),KY(end,:),'k'),plot(KX(:,1),KY(:,1),'k'),plot(KX(:,end-1),KY(:,end-1),'k')
text(0.5,0.9 ,'Re(B)','Units', 'normalized','HorizontalAlignment', 'center', 'FontSize', 10);
text(0.6,1.1 ,'$m<n, a<b$','Units', 'normalized','HorizontalAlignment', 'center', 'FontSize', 14,'Color','r',  'Interpreter', 'latex');
text(0.55,0.1 ,'k_1','Units', 'normalized', 'FontSize', 10);text(0.1,0.35 ,'k_2','Units', 'normalized', 'FontSize', 10);
text(KX(1,1)-0.3,KY(1,1)-0.3,'0', 'FontSize', 10,'HorizontalAlignment', 'center')
text(KX(1,end-1)-0.1,KY(1,end-1)-0.4,'1', 'FontSize', 10,'HorizontalAlignment', 'center')
text(KX(end-1,1)-0.3,KY(end-1,1)-0.3,'1', 'FontSize', 10,'HorizontalAlignment', 'center')
text(xte, 1.1, '(d)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');

h7=axes_list(8);axes(h7);   % 将其设为当前激活的子图
pcolor(KX(1:end-1,1:end-1),KY(1:end-1,1:end-1),imag(B)),hold on
patch([-1.1 NaN],[1.1 NaN],[-4 NaN],'EdgeColor','interp'),patch([-1.1 NaN],[1.1 NaN],[4 NaN],'EdgeColor','interp')
shading interp,xlabel('k_1'),ylabel('k_2'),axis equal,axis off,colormap(h7,b_w_r)
plot(KX(1,:),KY(1,:),'k'),plot(KX(end,:),KY(end,:),'k'),plot(KX(:,1),KY(:,1),'k'),plot(KX(:,end-1),KY(:,end-1),'k')
text(0.5,0.9 ,'Im(B)','Units', 'normalized','HorizontalAlignment', 'center', 'FontSize', 10);
text(0.55,0.1 ,'k_1','Units', 'normalized', 'FontSize', 10);text(0.1,0.35 ,'k_2','Units', 'normalized', 'FontSize', 10);
text(KX(1,1)-0.3,KY(1,1)-0.3,'0', 'FontSize', 10,'HorizontalAlignment', 'center')
text(KX(1,end-1)-0.1,KY(1,end-1)-0.4,'1', 'FontSize', 10,'HorizontalAlignment', 'center')
text(KX(end-1,1)-0.3,KY(end-1,1)-0.3,'1', 'FontSize', 10,'HorizontalAlignment', 'center')

t1=1;a=0.2;c=0.1;b=0.1;d=0.2;
[C,k11,w,B] = Chern_mol(t1,a,b,c,d,L);C
h6=axes_list(4);axes(h6);   % 将其设为当前激活的子图
pcolor(KX(1:end-1,1:end-1),KY(1:end-1,1:end-1),real(B)),hold on
patch([-1.1 NaN],[1.1 NaN],[-4 NaN],'EdgeColor','interp'),patch([-1.1 NaN],[1.1 NaN],[4 NaN],'EdgeColor','interp')
shading interp,xlabel('k_1'),ylabel('k_2'),axis equal,axis off,colormap(h6,b_w_r)
plot(KX(1,:),KY(1,:),'k'),plot(KX(end,:),KY(end,:),'k'),plot(KX(:,1),KY(:,1),'k'),plot(KX(:,end-1),KY(:,end-1),'k')
text(0.5,0.9 ,'Re(B)','Units', 'normalized','HorizontalAlignment', 'center', 'FontSize', 10);
text(0.6,1.1 ,'$m>n, a<b$','Units', 'normalized','HorizontalAlignment', 'center', 'FontSize', 14,'Color','r',  'Interpreter', 'latex');
text(0.55,0.1 ,'k_1','Units', 'normalized', 'FontSize', 10);text(0.1,0.35 ,'k_2','Units', 'normalized', 'FontSize', 10);
text(KX(1,1)-0.3,KY(1,1)-0.3,'0', 'FontSize', 10,'HorizontalAlignment', 'center')
text(KX(1,end-1)-0.1,KY(1,end-1)-0.4,'1', 'FontSize', 10,'HorizontalAlignment', 'center')
text(KX(end-1,1)-0.3,KY(end-1,1)-0.3,'1', 'FontSize', 10,'HorizontalAlignment', 'center')
text(xte, 1.1, '(e)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');

h7=axes_list(9);axes(h7);   % 将其设为当前激活的子图
pcolor(KX(1:end-1,1:end-1),KY(1:end-1,1:end-1),imag(B)),hold on
patch([-1.1 NaN],[1.1 NaN],[-4 NaN],'EdgeColor','interp'),patch([-1.1 NaN],[1.1 NaN],[4 NaN],'EdgeColor','interp')
shading interp,xlabel('k_1'),ylabel('k_2'),axis equal,axis off,colormap(h7,b_w_r)
plot(KX(1,:),KY(1,:),'k'),plot(KX(end,:),KY(end,:),'k'),plot(KX(:,1),KY(:,1),'k'),plot(KX(:,end-1),KY(:,end-1),'k')
text(0.5,0.9 ,'Im(B)','Units', 'normalized','HorizontalAlignment', 'center', 'FontSize', 10);
text(0.55,0.1 ,'k_1','Units', 'normalized', 'FontSize', 10);text(0.1,0.35 ,'k_2','Units', 'normalized', 'FontSize', 10);
text(KX(1,1)-0.3,KY(1,1)-0.3,'0', 'FontSize', 10,'HorizontalAlignment', 'center')
text(KX(1,end-1)-0.1,KY(1,end-1)-0.4,'1', 'FontSize', 10,'HorizontalAlignment', 'center')
text(KX(end-1,1)-0.3,KY(end-1,1)-0.3,'1', 'FontSize', 10,'HorizontalAlignment', 'center')

t1=1;a=0.1;c=0.1;b=0.25;d=0.1;
[C,k11,w,B] = Chern_mol(t1,a,b,c,d,L);C
h8=axes_list(5);axes(h8);   % 将其设为当前激活的子图
pcolor(KX(1:end-1,1:end-1),KY(1:end-1,1:end-1),real(B)),hold on
patch([-1.1 NaN],[1.1 NaN],[-4 NaN],'EdgeColor','interp'),patch([-1.1 NaN],[1.1 NaN],[4 NaN],'EdgeColor','interp')
shading interp,xlabel('k_1'),ylabel('k_2'),axis equal,axis off,colormap(h8,b_w_r)
plot(KX(1,:),KY(1,:),'k'),plot(KX(end,:),KY(end,:),'k'),plot(KX(:,1),KY(:,1),'k'),plot(KX(:,end-1),KY(:,end-1),'k')
text(0.5,0.9 ,'Re(B)','Units', 'normalized','HorizontalAlignment', 'center', 'FontSize', 10);
text(0.6,1.1 ,'$m=n, a>b$','Units', 'normalized','HorizontalAlignment', 'center', 'FontSize', 14,'Color','r',  'Interpreter', 'latex');
text(0.55,0.1 ,'k_1','Units', 'normalized', 'FontSize', 10);text(0.1,0.35 ,'k_2','Units', 'normalized', 'FontSize', 10);
text(KX(1,1)-0.3,KY(1,1)-0.3,'0', 'FontSize', 10,'HorizontalAlignment', 'center')
text(KX(1,end-1)-0.1,KY(1,end-1)-0.4,'1', 'FontSize', 10,'HorizontalAlignment', 'center')
text(KX(end-1,1)-0.3,KY(end-1,1)-0.3,'1', 'FontSize', 10,'HorizontalAlignment', 'center')
text(xte, 1.1, '(f)', 'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');

h9=axes_list(10);axes(h9);   % 将其设为当前激活的子图
pcolor(KX(1:end-1,1:end-1),KY(1:end-1,1:end-1),imag(B)),hold on
patch([-1.1 NaN],[1.1 NaN],[-4 NaN],'EdgeColor','interp'),patch([-1.1 NaN],[1.1 NaN],[4 NaN],'EdgeColor','interp')
shading interp,xlabel('k_1'),ylabel('k_2'),axis equal,axis off,colormap(h9,b_w_r)
plot(KX(1,:),KY(1,:),'k'),plot(KX(end,:),KY(end,:),'k'),plot(KX(:,1),KY(:,1),'k'),plot(KX(:,end-1),KY(:,end-1),'k')
text(0.5,0.9 ,'Im(B)','Units', 'normalized','HorizontalAlignment', 'center', 'FontSize', 10);
text(0.55,0.1 ,'k_1','Units', 'normalized', 'FontSize', 10);text(0.1,0.35 ,'k_2','Units', 'normalized', 'FontSize', 10);
text(KX(1,1)-0.3,KY(1,1)-0.3,'0', 'FontSize', 10,'HorizontalAlignment', 'center')
text(KX(1,end-1)-0.1,KY(1,end-1)-0.4,'1', 'FontSize', 10,'HorizontalAlignment', 'center')
text(KX(end-1,1)-0.3,KY(end-1,1)-0.3,'1', 'FontSize', 10,'HorizontalAlignment', 'center')
c=colorbar;%set(c, 'Orientation', 'horizontal'); 
c.Position=[0.92,0.2,0.01,0.6];c.FontSize=12; % Position and size of the colorbar