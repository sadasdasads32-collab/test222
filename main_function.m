clear; clc; close all;
%% 双层 QZS 隔振系统：FRF + Level-1(定频扫力) + 稳定性可视化
% -------------------------------------------------------------------------
% 更新日志：
% 1. 幅值定义统一：强制使用基波幅值 A = sqrt(A1c^2 + B1s^2)
% 2. 扫力范围扩大：默认关注 Force 0 ~ 12 范围
% 3. 稳定性可视化：自动检测 S 形曲线的拐点 (Limit Points)，将中间段标为虚线
% -------------------------------------------------------------------------

global Fw FixedOmega
FixedOmega = [];

%% 0. 统一幅值定义 (Function Handle)
% 输入 x 为系统状态矩阵 (15 x N) 或向量
% 状态索引: 2->A1c (主质量余弦), 3->B1s (主质量正弦)
% 状态索引: 7->A2c (副质量余弦), 8->B2s (副质量正弦)
calc_amp_primary   = @(x) sqrt(x(2,:).^2 + x(3,:).^2);
calc_amp_secondary = @(x) sqrt(x(7,:).^2 + x(8,:).^2);

%% 1. 系统参数
mu     = 0.1;
al1    = 0.01;
ga1    = 100.0;     % 非线性强度
ze1    = 0.005;
be1    = 1.0;
be2    = 1.0;
lam    = 0.5;
sig    = 0.1;
kap_e  = 1.0;

Fw = 1.0; % 初始激励力
sysP = [be1, be2, mu, al1, ga1, ze1, lam, sig, kap_e];

%% 2. 计算 FRF (幅频响应)
fprintf('--------------------------------------------------\n');
fprintf('STEP 1: 计算基础幅频响应曲线 (FRF)...\n');
FixedOmega = []; % 确保是扫频模式
x_mrc = FRF(sysP); % 运行 FRF 求解

W_frf  = x_mrc(end,:);
A1_frf = calc_amp_primary(x_mrc);
A2_frf = calc_amp_secondary(x_mrc);

% --- 自动寻找 Level-1 的最佳起点 ---
% 策略：找到主共振峰值，并在其右侧（硬特性折叠区）选一个频率
[max_A, i_pk] = max(A1_frf);
W_pk = W_frf(i_pk);
target_W = W_pk + 0.5; % 稍微偏右，更容易扫出 S 形
[~, idx_L1] = min(abs(W_frf - target_W));
x_L1_start = x_mrc(:, idx_L1)'; % 提取该频率下的状态作为 L1 初值

fprintf('FRF 计算完成。峰值频率 W=%.3f, 峰值幅值 A=%.3f\n', W_pk, max_A);
fprintf('已自动选择 Level-1 扫描频率: Omega = %.3f\n', x_L1_start(end));

%% 3. 绘图 1: FRF 概览
figure(1); clf;
subplot(1,2,1);
plot(W_frf, A1_frf, 'b-', 'LineWidth', 1.5); grid on;
hold on; plot(x_L1_start(end), calc_amp_primary(x_L1_start'), 'ro', 'MarkerFaceColor','r');
title('Primary Mass FRF (Fundamental)'); xlabel('\Omega'); ylabel('Amp A_1');
text(x_L1_start(end), calc_amp_primary(x_L1_start'), '  \leftarrow L1 Start', 'FontSize', 8);

subplot(1,2,2);
plot(W_frf, A2_frf, 'k-', 'LineWidth', 1.5); grid on;
title('Secondary Mass FRF'); xlabel('\Omega'); ylabel('Amp A_2');

%% 4. 执行 Level-1 单频扫力 (含稳定性可视化)
% 这是一个演示“论文级”图表生成的过程
fprintf('--------------------------------------------------\n');
fprintf('STEP 2: 执行单频 Level-1 扫描 (Force Response)...\n');

do_single = input('是否开始计算单频 Force-Amplitude 曲线? (y/n, 推荐 y): ', 's');
if isempty(do_single), do_single = 'y'; end

if strcmpi(do_single, 'y')
    % 自定义步长和参数
    force_step = 0.005; 
    
    % 调用 L1 追踪
    xt = L1(x_L1_start, sysP, force_step);
    
    if ~isempty(xt)
        % 提取数据
        Force_vec = xt(end, :);     % 最后一维是 Force
        Amp1_vec  = calc_amp_primary(xt);
        
        % --- 稳定性/分支识别算法 ---
        % 逻辑：在 S 形曲线上，拐点 (Limit Points) 处的切线 dF/dA = 0
        % 也就是 Force 对弧长的导数变号的位置
        
        % 1. 计算 Force 的差分
        dF = diff(Force_vec);
        
        % 2. 寻找符号变化点 (即拐点索引)
        % 忽略开头和结尾的微小波动
        sign_change_idx = find(diff(sign(dF)) ~= 0);
        
        % 过滤掉可能是数值噪声的微小拐点
        valid_bif_idx = [];
        for k = 1:length(sign_change_idx)
            idx = sign_change_idx(k);
            % 简单的过滤：要求该段长度大于 10 个点才算有效分支
            if idx > 10 && idx < length(Force_vec)-10
                valid_bif_idx = [valid_bif_idx, idx];
            end
        end
        
        % 绘图准备
        figure(4); clf; hold on; grid on;
        title_str = sprintf('Force Response (\\Omega=%.2f) with Stability', x_L1_start(end));
        title(title_str); 
        xlabel('Excitation Force F'); ylabel('Amplitude A_1 (Fundamental)');
        
        if isempty(valid_bif_idx)
            % 如果没检测到拐点（单调曲线），直接画实线
            plot(Force_vec, Amp1_vec, 'b-', 'LineWidth', 2);
            legend('Stable Solution');
        else
            % 如果检测到拐点，分段画
            % 通常顺序是：稳定(下) -> 不稳定(中) -> 稳定(上)
            
            % 第一段：起点到第一个拐点 (通常是稳定)
            idx1 = valid_bif_idx(1);
            p1 = plot(Force_vec(1:idx1), Amp1_vec(1:idx1), 'b-', 'LineWidth', 2);
            
            % 中间段：拐点之间 (通常是不稳定)
            if length(valid_bif_idx) >= 2
                idx2 = valid_bif_idx(2);
                p2 = plot(Force_vec(idx1:idx2), Amp1_vec(idx1:idx2), 'r--', 'LineWidth', 2);
                
                % 第三段：第二个拐点之后 (通常是稳定)
                p3 = plot(Force_vec(idx2:end), Amp1_vec(idx2:end), 'b-', 'LineWidth', 2);
                
                % 标注拐点 (Limit Points)
                plot(Force_vec(idx1), Amp1_vec(idx1), 'ko', 'MarkerFaceColor', 'k');
                plot(Force_vec(idx2), Amp1_vec(idx2), 'ko', 'MarkerFaceColor', 'k');
                text(Force_vec(idx1), Amp1_vec(idx1), ' LP1', 'VerticalAlignment','bottom');
                text(Force_vec(idx2), Amp1_vec(idx2), ' LP2', 'VerticalAlignment','top');
                
                legend([p1, p2], {'Stable', 'Unstable'}, 'Location', 'best');
            else
                % 只有一个拐点的情况
                plot(Force_vec(idx1:end), Amp1_vec(idx1:end), 'r--', 'LineWidth', 2);
            end
        end
        
        xlim([0, max(Force_vec)*1.1]);
        fprintf('绘图完成。红色虚线表示潜在的不稳定分支 (Unstable Branch)。\n');
    end
end

%% 5. 多频率扫描 (Waterfall / Parametric Study)
fprintf('--------------------------------------------------\n');
fprintf('STEP 3: 多频率 Level-1 扫描 (Figure 3)...\n');
do_multi = input('是否更新多频率汇总图? (y/n): ', 's');

if strcmpi(do_multi, 'y')
    figure(3); clf; hold on; grid on;
    title('Force-Amplitude Curves at Different Frequencies');
    xlabel('Force F'); ylabel('Amplitude A_1');
    
    % 选取几个特征频率
    % 包含：线性区(远离峰值), 峰值左侧, 峰值处, 峰值右侧(强非线性)
    target_Omegas = [W_pk-2, W_pk-0.5, W_pk, W_pk+0.5, W_pk+2];
    
    % 颜色映射
    colors = lines(length(target_Omegas));
    
    legend_str = {};
    
    for i = 1:length(target_Omegas)
        Wi = target_Omegas(i);
        % 找最近的初值
        [~, id] = min(abs(W_frf - Wi));
        x_start = x_mrc(:, id)';
        
        fprintf('  正在扫描 Omega = %.2f ...\n', Wi);
        xt_i = L1(x_start, sysP, 0.005); % 步长 0.005
        
        if ~isempty(xt_i)
            F_i = xt_i(end,:);
            A_i = calc_amp_primary(xt_i);
            
            % 这里只画细实线，不标稳定性，避免图太乱
            plot(F_i, A_i, 'Color', colors(i,:), 'LineWidth', 1.5);
            legend_str{end+1} = sprintf('\\Omega=%.2f', Wi);
        end
    end
    legend(legend_str, 'Location', 'bestoutside');
end

fprintf('--------------------------------------------------\n');
fprintf('全部分析完成。\n');