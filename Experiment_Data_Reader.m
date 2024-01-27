clc, clear, close all;
P_=xlsread('D:\Users\Jianzhe Wei\项目研究\L3_0.18_R.xlsx');
P_=P_';
[num,txt,raw]=xlsread('D:\Users\Jianzhe Wei\项目研究\cells_trajectory_lambda=1.60_v2.xlsx');
cell_index = unique(raw(2:end, 1)); % 从第二行开始，取第一列的数据，然后去除重复的值
n = length(cell_index); % 计算有多少种 Cell_index，这里应该是 4
Q = cell(n, 1); % 创建一个 n*1 的 cell 数组，每个元素都是一个 cell
len = size(raw, 1);
for i = 2:len % 从第二行开始，遍历 data 的每一行
    ci = raw{i, 1}; % 取出当前行的 Cell_index 的值
    t = raw{i, 2}; % 取出当前行的 Time 的值
    g = raw{i, 3}; % 取出当前行的 green 的值
    r = raw{i, 4}; % 取出当前行的 Red 的值
    j = find(strcmp(cell_index, ci)); % 在 cell_index 中找到当前行的 Cell_index 对应的位置，这个位置就是 Q 的第一维的索引
    Q{j} = [Q{j}; t, g, r]; % 将 Time，green，Red 的数据追加到 Q 的第 j 个 cell 中
end



%%
% clc, clear, close all;
% Q=importdata('D:\Program Files\Polyspace\R2021a\bin\Research_Projects\Data_Experiment\P([t,G,R])_trajectories_lambda=1.60_Experimental.mat');
% n = length(Q); % 计算 Q 的长度，这里应该是 4
% mask = false(n, 1); % 创建一个和 Q 的长度相同的逻辑数组，用来标记哪些 cell 需要被删除
% for i = 1:n % 遍历 Q 中的每个 cell
%     if size(Q{i}, 1) ~= 90 % 如果当前 cell 的行数不等于 2，说明需要删除
%         mask(i) = true; % 将当前 cell 的 mask 设为 true
%     end
% end
% Q(mask) = []; % 删除 mask 为 true 的 cell
% 
% for k=1:1358
%     P_(k,:,:)=Q{k}';
% end
% t0_(:,:)=P_(:,1,:); t0_=t0_-t0_(:,1);
% x0_=P_(:,2:3,:);
% for j=1:1358
%     t0(:)=t0_(j,:); x0(:,:)=x0_(j,:,:);
%     x=observed_time_trans(t0,x0,13.9,0.1);
%     P(j,:,:)=[x(2,:);x(1,:)];
% end
% 
% function x=observed_time_trans(t0,x0,T,dt)
% i=1; j=1; i0=T/dt+1; t=0:dt:T;
% x=zeros(2,i0);
% while (i <= i0)
%     if t(i) >= t0(j)
%         j=j+1;
%         continue
%     elseif (t(i) < t0(j))
%         x(:,i)=x0(:,j-1);
%         i=i+1;
%     end
% end
% end





%%
clc, clear, close all;
P=importdata('D:\Program Files\Polyspace\R2021a\bin\Research_Projects\Data_Experiment\P([R,G])_SteadyState_lambda=0.18_Experimental.mat');
% P=log2(P_);
% P=0.7.*P_;
plot(P(1,:),P(2,:),'o');
xlabel('$[R]$','Interpreter','latex')
ylabel('$[G]$','Interpreter','latex','Rotation',0)
title('Scatter in steady state, lambda=0.18, Experimental');

figure(2)
[x,y,h]=heatmap(P,-30,12030,60,-0.4,70.4,0.8);
pcolor(x,y,h);
colormap('turbo');
colorbar;
shading interp;
axis on;
xlabel('$[R]$','Interpreter','latex')
ylabel('$[G]$','Interpreter','latex','Rotation',0)
title('Heatmap in steady state, lambda=0.18, Experimental');

R=P(1,:).*317/1167; G=P(2,:).*11/12.2;
[xR,hR]=distribution(R,0,3500,10); Distribution_R=[xR',hR'];
[xG,hG]=distribution(G,0,70,1);  Distribution_G=[xG',hG'];
figure(3)
% semilogx(xR,hR,'linewidth',1.5);
plot(xR,hR,'linewidth',1.5);
title('Distribution of R in steady state, lambda=0.18, Experimental');
figure(4)
% semilogx(xG,hG,'linewidth',1.5);
plot(xG,hG,'linewidth',1.5);
title('Distribution of G in steady state, lambda=0.18, Experimental');
mean_R=mean(R);
mean_G=mean(G);
var_R=var(R);
var_G=var(G);
CV_R=sqrt(var_R)/mean_R; 
CV_G=sqrt(var_G)/mean_G;

function [x,h]=distribution(P,Xl,Xr,dx)
x=Xl+dx/2:dx:Xr-dx/2;
edges = Xl:dx:Xr;
h=histcounts(P,edges,'Normalization','pdf');
end

function [x_,y_,h_]=heatmap(P,Xl,Xr,dx,Yl,Yr,dy)
Xedges = Xl:dx:Xr;
Yedges = Yl:dy:Yr;
h=histcounts2(P(1,:),P(2,:),Xedges,Yedges,'Normalization','pdf');
h_=h';
x=Xl+dx/2:dx:Xr-dx/2;
y=Yl+dy/2:dy:Yr-dy/2;
[x_,y_]=meshgrid(x,y);
end