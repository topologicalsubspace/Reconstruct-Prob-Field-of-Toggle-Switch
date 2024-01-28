clc, clear, close all;
%% ------实验数据预处理，得到结果为P_new_-------------------------------------
%%
Q=importdata('D:\Program Files\Polyspace\R2021a\bin\Research_Projects\Data_Experiment\P([t,G,R])_trajectories_lambda=1.60_Experimental.mat');
n = length(Q); 
mask = false(n, 1); 
for i = 1:n 
    if size(Q{i}, 1) ~= 90 
        mask(i) = true; 
    end
end
Q(mask) = []; 
for k=1:1358
    P(k,:,:)=Q{k}';
end

% 看实验拍摄的∆𝑡分布(图)
t(:,:)=P(:,1,:);
difft = diff(t,1,2);
diff_t = difft(:);
[x_R,h_R]=distribution(diff_t,0.15,0.165,0.0001);
plot(x_R,h_R,'o','linewidth',1.5);
xlabel('$\Delta t$','Interpreter','latex')
ylabel('$Probability(Normalized)$','Interpreter','latex')
title('Probability Density Function of \Deltat');

% 根据∆𝑡分布(图)合理设置t的间隔dt，重新得到P_new_
global dt
dt = 0.313; % 你需要设置 dt 的值
x0_=P(:,2:3,:);
for j=1:1358
    t0(:)=t(j,:); x0(:,:)=x0_(j,:,:);
    [t_new,x_R]=observed_time_trans(t0,x0,13.9,dt);
    P_new(j,:,:)=[x_R(2,:);x_R(1,:)];
end
P_new_ = squeeze(P_new(:,:,2:end));

%% ------计算并估计合适间隔 ∆𝑡, ∆𝑅, ∆𝐺 -------------------------------------
%%
R(:,:)=P_new_(:,1,:); G(:,:)=P_new_(:,2,:);
diffR = diff(R,1,2); diffG = diff(G,1,2);
%% 画总体分布图
figure(1);
diff_R = diffR(:);
[x_R,h_R]=distribution(diff_R,-50,50,0.1);
plot(x_R,h_R,'linewidth',2);
xticks(-50:5:50); % 设置 x 轴的刻度
xlabel('$\Delta R$','Interpreter','latex')
ylabel('$Probability (Normalized)$','Interpreter','latex')
title('Probability Density Function of \DeltaR');

figure(2);
diff_G = diffG(:);
[x_G,h_G]=distribution(diff_G,-10,10,0.05);
plot(x_G,h_G,'linewidth',2);
xticks(-10:1:10); % 设置 x 轴的刻度
xlabel('$\Delta G$','Interpreter','latex')
ylabel('$Probability (Normalized)$','Interpreter','latex')
title('Probability Density Function of \DeltaG');
%% 画散点图
[n, ~, len] = size(P_new_);
figure(3); 
hold on; 
for i = 1:len-1
    R_t=[i*dt*ones(n,1),diffR(:,i)];
    scatter(R_t(:,1),R_t(:,2)); 
    mediansR(i) = median(diffR(:,i)); 
end
pR = plot(dt*(1:len-1), mediansR,'r','linewidth',2); 
legend(pR, 'Median'); 
hold off; 
axis([0,14,-120,120]);
xlabel('$t$','Interpreter','latex')
ylabel('$\Delta R$','Interpreter','latex')
title('Scatter of \DeltaR');

figure(4); 
hold on; 
for i = 1:len-1
    G_t=[i*dt*ones(n,1),diffG(:,i)];
    scatter(G_t(:,1),G_t(:,2)); 
    mediansG(i) = median(diffG(:,i)); 
end
pG = plot(dt*(1:len-1), mediansG,'r','linewidth',2);
legend(pG, 'Median'); 
hold off; 
axis([0,14,-50,50]);
xlabel('$t$','Interpreter','latex')
ylabel('$\Delta G$','Interpreter','latex')
title('Scatter of \DeltaG');

%% ------正式构建场-------------------------------------------------------------
%%

% 接下来考虑构建场，见Part " Reconstruct_Prob_Field_2_main "

%%
%%
function [t,x]=observed_time_trans(t0,x0,T,dt) 
    i=1; j=1; i0=floor(T/dt+1); t=0:dt:T; 
    x=zeros(2,i0); 
    while (i <= i0) 
        if t(i) < t0(1) 
            x(:,i) = NaN; % 当t(i) < t0(1)时，将x(:,i)设置为NaN
            i = i + 1;
        elseif t(i) >= t0(j) 
            j=j+1; 
        elseif (t(i) < t0(j)) 
            x(:,i)=x0(:,j-1); 
            i=i+1; 
        end 
    end 
end

function [x,h]=distribution(P,Xl,Xr,dx)
x=Xl+dx/2:dx:Xr-dx/2;
edges = Xl:dx:Xr;
h=histcounts(P,edges,'Normalization','pdf');
end