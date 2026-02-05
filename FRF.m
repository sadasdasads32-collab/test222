function [x] = FRF(sysP)
%% 使用谐波平衡法 (HBM) 和弧长连续法追踪幅频响应曲线 (FRF)
% 适配 2-DOF QZS 系统 (15 维变量)
global Fw FixedOmega
FixedOmega = []; % 确保进入 FRF 模式（mu 代表频率）

% 1. 初始化 15 维解向量
% y = [xi1(A0, A1c, B1s, A3c, B3s), xi2(A0, A1c, B1s, A3c, B3s), Q(A0, A1c, B1s, A3c, B3s)]
y = zeros(15, 1); 
% 给出一个微小的初值防止雅可比矩阵在纯零点奇异
y(2) = 0.1; % 给 xi1 一个初始的余弦项幅值

tic
fprintf('开始执行 FRF 扫频...\n');

% 2. 确定追踪步数 (Branch)
% 根据激励力 Fw 动态调整，Fw 越大非线性越强，需要的步数越多
if Fw < 0.2
    branch = 2000;
elseif Fw < 1.1
    branch = 5000;
else
    branch = 10000;
end

% 3. 设定起始频率 (无量纲频率 Omega)
% 从低频 0.2 开始，跨越到 0.25 获取初始切线
omega0 = 0.2;
omega1 = 0.21;

% 初始点计算：寻找频率为 omega0 时的稳态解（Newton 固定 omega）
y_init = [y; omega0];
x0_full = newton('nondim_temp2', y_init, sysP);

% 第二点计算：寻找频率为 omega1 时的稳态解，为弧长连续提供初始切线
y_init1 = [x0_full(1:15); omega1];
x1_full = newton('nondim_temp2', y_init1, sysP);

% 重要：弧长连续的 tracking_file（nondim_temp2）在 continuation 模式下只接收 15维状态
x0 = x0_full(1:15);
x1 = x1_full(1:15);

% 4. 执行弧长连续法追踪路径
% 'nondim_temp2' 会返回残差、雅可比方向及约束
% 使用自动化版弧长连续（不包含旧版的 temp/Floquet 依赖）
[x, ~] = branch_follow2('nondim_temp2', branch, omega0, omega1, x0, x1, sysP);

% 5. 检查是否追踪到了预定频率范围 (例如 Omega = 3.5)
last_omega = x(end, end);
fprintf('FRF 追踪完成。最终频率: %.2f\n', last_omega);

toc
end