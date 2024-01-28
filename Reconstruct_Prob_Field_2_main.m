clc, clear, close all;
P_new_=importdata('P_OU_L_dt=0.1.mat');
%% ------计算 Prob. rate, 即transition_rates_xr_cells_new-------------------
%%
global dt
dt = 0.1; % 设置 dt 的值
[n, ~, len] = size(P_new_);
grid_size = [4, 4];  % 设置网格大小
bundary_x_l = 100; bundary_x_r = 200;  % 设置 x 轴边界（统计全局）
bundary_y_l = 100; bundary_y_r = 200;  % 设置 y 轴边界（统计全局） 
num_grids = [(bundary_x_r - bundary_x_l)/grid_size(1), (bundary_y_r - bundary_y_l)/grid_size(2)];  
max_tau = len - 1;  % 设置最大的时间间隔
threshold = 3;  % 设置阈值
transition_rates_xr = zeros(num_grids(1), num_grids(2), len-1, max_tau);

delete(gcp('nocreate'))
parpool('local',14)
parfor tau = 1:max_tau
    transition_counts = zeros(num_grids(1), num_grids(2), len-1);
    grid_counts = zeros(num_grids(1), num_grids(2), len);
    for t = 1:len-tau
        for i = 1:n
            x = P_new_(i, 1, t);
            y = P_new_(i, 2, t);
            next_x = P_new_(i, 1, t+tau);
            next_y = P_new_(i, 2, t+tau);
            
            % 如果 (x, y) 或 (next_x, next_y) 超过全局边界，跳过当前的样本
            if x < bundary_x_l || x > bundary_x_r || y < bundary_y_l || y > bundary_y_r || next_x < bundary_x_l || next_x > bundary_x_r || next_y < bundary_y_l || next_y > bundary_y_r
                continue;
            end
            
            grid_x = min(floor((x - bundary_x_l)/grid_size(1)) + 1, num_grids(1));
            grid_y = min(floor((y - bundary_y_l)/grid_size(2)) + 1, num_grids(2));
            next_grid_x = min(floor((next_x - bundary_x_l)/grid_size(1)) + 1, num_grids(1));
            next_grid_y = min(floor((next_y - bundary_y_l)/grid_size(2)) + 1, num_grids(2));
            
            % 计算每个网格中的样本数目
            grid_counts(grid_x, grid_y, t) = grid_counts(grid_x, grid_y, t) + 1;
            
            % 检查在 t 和 t+tau 之间的所有时间步，样本是否都在原网格或相邻网格内
            in_between = true;
            for t_prime = t+1:t+tau-1
                x_prime = P_new_(i, 1, t_prime);
                y_prime = P_new_(i, 2, t_prime);
                grid_x_prime = min(floor((x_prime - bundary_x_l)/grid_size(1)) + 1, num_grids(1));
                grid_y_prime = min(floor((y_prime - bundary_y_l)/grid_size(2)) + 1, num_grids(2));
                % 注意if中设置转移的方向
%                 if ~((grid_x_prime == grid_x || grid_x_prime == grid_x + 1) && grid_y_prime == grid_y)   % x轴正向
%                 if ~((grid_x_prime == grid_x || grid_x_prime == grid_x - 1) && grid_y_prime == grid_y)   % x轴反向
%                 if ~((grid_x_prime == grid_x) && (grid_y_prime == grid_y || grid_y_prime == grid_y + 1))   % y轴正向
                if ~((grid_x_prime == grid_x) && (grid_y_prime == grid_y || grid_y_prime == grid_y - 1))   % y轴反向
                    
                    in_between = false;
                    break;
                end
            end
            
            % 如果在 t 和 t+tau 之间的所有时间步，样本都在原网格或相邻网格内，那么计算转移计数
            % 注意if中设置转移的方向
%             if in_between && next_grid_x == grid_x + 1 && next_grid_y == grid_y   % x轴正向
%             if in_between && next_grid_x == grid_x - 1 && next_grid_y == grid_y   % x轴反向
%             if in_between && next_grid_x == grid_x && next_grid_y == grid_y + 1   % y轴正向
            if in_between && next_grid_x == grid_x && next_grid_y == grid_y - 1   % y轴反向

                transition_counts(grid_x, grid_y, t) = transition_counts(grid_x, grid_y, t) + 1;
            end
        end    
    end

    % 计算转移速率
    transition_rates = transition_counts ./ grid_counts(:, :, 1:end-1);

    % 如果原本的网格格子中样本数目（即基数）小于阈值，将转移速率设置为 nan
    transition_rates(grid_counts(:, :, 1:end-1) < threshold) = NaN;
    transition_rates_xr(:,:,:,tau) = transition_rates;
end
delete(gcp('nocreate'))

threshold2 = 1;  % 阈值
% 计算每个网格中非 NaN 的个数
num_non_nan = sum(~isnan(transition_rates_xr), [3, 4]);

transition_rates_xr_cells = cell(size(transition_rates_xr, 1), size(transition_rates_xr, 2));
for i = 1:size(transition_rates_xr, 1)
    for j = 1:size(transition_rates_xr, 2)
        % 如果 NaN 的个数大于等于阈值，那么将 transition_rates_xr(i, j, :, :) 存储到 cell 数组中
        if num_non_nan(i, j) > threshold2
            transition_rates_xr_cells{i, j} = squeeze(transition_rates_xr(i, j, :, :));
        end
    end
end

% 删除 transition_rates_xr, 减少运行内存开支
clear transition_rates_xr

notEmptyCells = ~cellfun(@isempty, transition_rates_xr_cells);
notEmptyCellsIndices = find(notEmptyCells)';
notEmptyCellsElements = transition_rates_xr_cells(notEmptyCellsIndices);
new_transition_rates_xr_cells = cell(size(notEmptyCellsElements));
transition_rates_xr_cells_new = transition_rates_xr_cells;

% 对每个非空 cell 进行操作
delete(gcp('nocreate'))
parpool('local',14)
parfor idx = 1:length(notEmptyCellsIndices)
    A = notEmptyCellsElements{idx};    
    B = [];    
    for tau = 1:size(A, 2)
        % 获取当前列的有效部分（去除 nan）
        column = A(:, tau);
        validColumn = column(~isnan(column));       
        % 将 tau*dt 与当前列的有效部分组合成一个新的数组，并添加到新的二维矩阵中
        B = [B; [repmat(tau*dt, length(validColumn), 1), validColumn]];        
    end    
    new_transition_rates_xr_cells{idx} = B;
end
delete(gcp('nocreate'))

transition_rates_xr_cells_new(notEmptyCellsIndices) = new_transition_rates_xr_cells;

%% ------对初步Prob. rate进行分析--------------------------------------------
%% 观察𝑃(𝜏)的分布有多符合指数分布，通过fit_results_cells[R^2,拟合参数]判定
notEmptyCells = ~cellfun(@isempty, transition_rates_xr_cells_new);
notEmptyCellsIndices = find(notEmptyCells)';
notEmptyCellsElements = transition_rates_xr_cells_new(notEmptyCellsIndices);
new_fit_results_cells = cell(size(notEmptyCellsElements));
fit_results_cells = transition_rates_xr_cells_new;

% 定义拟合函数
f = fittype('exp(-b*(x-h))', 'independent', 'x', 'dependent', 'y');


delete(gcp('nocreate'))
parpool('local',14)
parfor idx = 1:length(notEmptyCellsIndices)
    B = notEmptyCellsElements{idx};
    % 进行拟合
    [fitobject, gof] = fit(B(:, 1), B(:, 2), f);
    % 将 R^2 值和拟合参数存储到新的 cell 中
    new_fit_results_cells{idx} = [gof.rsquare, fitobject.b, fitobject.h];
end
delete(gcp('nocreate'))

fit_results_cells(notEmptyCellsIndices) = new_fit_results_cells;

%% 做𝑃(𝜏)的分布(图)，以分析和确定P_ij的取值
i=5; j=13;
% b=(1/150)/(2*20); h=0;

B = transition_rates_xr_cells_new{i, j};

% 绘制散点图
scatter(B(:, 1), B(:, 2));
hold on;
% B_ss=exp(-b.*(B(:, 1)-h));
% plot(B(:, 1), B_ss,'r','linewidth',2);
% axis([0,14,0,1]);

% 绘制中位数线
unique_x = unique(B(:, 1));
median_y = arrayfun(@(x) median(B(B(:, 1) == x, 2)), unique_x);
plot(unique_x, median_y, 'r', 'LineWidth', 2);

xlabel('\tau');
ylabel('transition Prob.');
title(sprintf('transition Prob. for lattice (%d, %d)', i, j));
% 设置 y 轴的范围
ylim([0 1]);

%% ------得到P_ij的取值，即transition_rates_xr_new---------------------------
%%
[s1, s2] = size(transition_rates_xr_cells_new);
transition_rates_xr_new = nan(s1, s2);

delete(gcp('nocreate'))
parpool('local',14)
parfor i = 1:s1
    for j = 1:s2
        C = transition_rates_xr_cells_new{i, j};
        % 如果 B 不为空，且 B(:, 1) 包含 dt，计算 dt 对应的 B(:, 2) 的中位数
        if ~isempty(C) && any(C(:, 1) == dt)
            dt_median = median(C(C(:, 1) == dt, 2));
            transition_rates_xr_new(i, j) = dt_median;  % 直接赋值为矩阵元素
        end
    end
end
delete(gcp('nocreate'))

%% ------将4个方向transition_rates_xr_new重组，得到重构的向量平均场------------
%%
% 接下来将4个方向transition_rates_xr_new重组，得到重构的向量平均场，见Part " Reconstruct_Prob_Field_3_result "