clc, clear, close all;
P_new_=importdata('P_new_dt=0.1565.mat');
%% ------计算转移矩阵M(t0+tau), 即transition_rates---------------------------
%%
global dt
dt = 0.1565; % 设置 dt 的值
tau=2;
[n, ~, len] = size(P_new_);
grid_size = [10, 1];  % 设置网格大小
bundary_x_l = 0; bundary_x_r = 240;  % 设置 x 轴边界（统计全局）
bundary_y_l = 0; bundary_y_r = 40;  % 设置 y 轴边界（统计全局） 
num_grids = [(bundary_x_r - bundary_x_l)/grid_size(1), (bundary_y_r - bundary_y_l)/grid_size(2)];  
num_grids_total = num_grids(1)*num_grids(2);
transition_rates = zeros(num_grids_total, num_grids_total, len-tau);
transition_counts = zeros(num_grids_total, num_grids_total);
grid_counts = zeros(num_grids_total, len-tau);
threshold = 10;  % 设置阈值

for t = 1:len-tau
    local_grid_counts = zeros(num_grids_total,1);
    local_transition_counts = zeros(num_grids_total, num_grids_total);
    for i = 1:n
        x = P_new_(i, 1, t);
        y = P_new_(i, 2, t);
        next_x = P_new_(i, 1, t+tau);
        next_y = P_new_(i, 2, t+tau);

        if x < bundary_x_l || x > bundary_x_r || y < bundary_y_l || y > bundary_y_r || next_x < bundary_x_l || next_x > bundary_x_r || next_y < bundary_y_l || next_y > bundary_y_r
            continue;
        end

        grid_x = min(floor((x - bundary_x_l)/grid_size(1)) + 1, num_grids(1));
        grid_y = min(floor((y - bundary_y_l)/grid_size(2)) + 1, num_grids(2));
        next_grid_x = min(floor((next_x - bundary_x_l)/grid_size(1)) + 1, num_grids(1));
        next_grid_y = min(floor((next_y - bundary_y_l)/grid_size(2)) + 1, num_grids(2));
        
        % 把平面网格展开成列向量
        grid_index = (grid_y-1)*num_grids(1) + grid_x;
        next_grid_index = (next_grid_y-1)*num_grids(1) + next_grid_x;
        
        local_grid_counts(grid_index) = local_grid_counts(grid_index) + 1;
        local_transition_counts(next_grid_index, grid_index) = local_transition_counts(next_grid_index, grid_index) + 1;
        
    end
    
    grid_counts(:,t) = local_grid_counts;
    transition_counts(:,:,t) = local_transition_counts;
end
for t = 1:len-tau
    for grid_index = 1:num_grids_total
        transition_rates(:, grid_index, t) = transition_counts(:, grid_index, t) ./ grid_counts(grid_index, t);
    end
end
transition_rates(repmat(grid_counts < threshold, [1, num_grids_total, 1])) = NaN;
transition_rates(isnan(transition_rates)) = 0;

%% ------回溯验证，利用 C-K Eq.----------------------------------------------
%%
figure(1);
for t = 1:len-tau-1
    transition_rates_t0 = transition_rates(:,:,t);
    grid_counts_t0 = grid_counts(:,t);
    grid_counts_t0_tau = transition_rates_t0 * grid_counts_t0;

    subplot(len-tau-1, 1, t);
    plot(1:num_grids_total, grid_counts(:,t+tau), 'r');
    hold on;
    plot(1:num_grids_total, grid_counts_t0_tau, 'b');
    legend('Real', 'M \cdot P');
    title(['t0 = ', num2str((t+1)*dt), 'h']);
    xlabel('Grid Index');
    ylabel('Count');
end
%%
% 某一个时刻拿出来单独画出来看一下
figure(2);
    t = 80;
    transition_rates_t0 = transition_rates(:,:,t);
    grid_counts_t0 = grid_counts(:,t);
    grid_counts_t0_tau = transition_rates_t0 * grid_counts_t0;

    subplot(3, 1, 3);
    plot(1:num_grids_total, grid_counts(:,t+tau), 'r-', 'LineWidth', 2);
    hold on;
    plot(1:num_grids_total, grid_counts_t0_tau, 'b--', 'LineWidth', 2);
    yticks(0:10:60); 
    legend('Real', 'M \cdot P');
    title(['t0 = ', num2str((t+1)*dt), 'h']);
    xlabel('Grid Index');
    ylabel('Count');
    
%% ------验证马尔可夫性，利用 M(𝑡_0+2𝑑𝑡)= [M(𝑡_0+𝑑𝑡)]^2.----------------------
%%
clc, clear, close all;
transition_rates_dt=importdata('transition_rates_t0+dt.mat');
transition_rates_2dt=importdata('transition_rates_t0+2dt.mat');

len_dt = size(transition_rates_dt, 3);
len_2dt = size(transition_rates_2dt, 3);
transition_rates_dt_squared = zeros(size(transition_rates_2dt));

% 对每个时间步进行计算
for t = 1:len_2dt
    % 计算 transition_rates_dt 的平方
    transition_rates_dt_squared(:,:,t) = transition_rates_dt(:,:,t) ^ 2;
end

%%
t0 = 1;
transition_rates_dt_squared_t0 = transition_rates_dt_squared(:,:,t0);
transition_rates_2dt_t0 = transition_rates_2dt(:,:,t0);

% 创建一个新的图形窗口
figure;

% 找到两个矩阵中的最小值和最大值
minValue = min(min(transition_rates_dt_squared_t0(:)), min(transition_rates_2dt_t0(:)));
maxValue = max(max(transition_rates_dt_squared_t0(:)), max(transition_rates_2dt_t0(:)));

% 在第一个子图中显示 transition_rates_dt_squared_t0
subplot(1, 2, 1);
imagesc(transition_rates_dt_squared_t0, [minValue maxValue]);
colorbar;
title('transition\_rates\_dt\_squared\_t0');

% 在第二个子图中显示 transition_rates_2dt_t0
subplot(1, 2, 2);
imagesc(transition_rates_2dt_t0, [minValue maxValue]);
colorbar;
title('transition\_rates\_2dt\_t0');