function x = L1(x_interest, sysP, varargin)
    %% Level 1 Continuation: 力-响应曲线追踪 (S形曲线)
    global Fw FixedOmega 
    
    % --- 1. 参数处理 ---
    current_F = Fw;
    step = 0.001; 
    
    if nargin >= 3
        input_val = varargin{1};
        if abs(input_val) > 0.1
            current_F = input_val; % 如果输入值较大，视为初始力
        else
            step = input_val;      % 如果输入值较小，视为步长
        end
    end
    
    % --- 2. 变量提取 ---
    % x_interest: 来自 FRF 的某一列/点，通常是 [15个傅里叶系数; Omega]
    % 这里取前 15 个系数作为 Level-1 的未知量
    x_coeffs = x_interest(1:15)';
    
    % 动态获取频率
    if length(x_interest) >= 31
        fixed_omega = x_interest(31); 
    else
        fixed_omega = x_interest(end); 
    end
    
    % 锁定频率，激活 nondim_temp2 的模式 B
    FixedOmega = fixed_omega; 
    
    % --- 3. 构造启动向量 (16维: 15系数 + 1力) ---
    % 注意：我们不再需要计算切向量 phi0，branch_follow 会自动处理
    y_start = [x_coeffs; current_F];
    
    % --- 4. 初始点修正 ---
    % newton 看到 16 维输入，FixedOmega 非空 -> 知道最后一位是参数(力)
    % 它会固定力，求解前 15 个系数的平衡
    try
    x0 = newton('nondim_temp2', y_start, sysP);
    catch
        fprintf('L1: 初始点修正失败，尝试直接启动...\n');
        x0 = y_start;
    end
    
    % --- 5. 执行追踪 ---
    branch_len = 2000; 
    mu0 = x0(end);      % 获取修正后的力
    mu1 = mu0 + step;   % 作为第二点的参数
    
    fprintf('L1: 开始追踪 (Freq=%.2f, Start Force=%.2f)...\n', fixed_omega, mu0);
    
    % branch_follow2 现在处理的是 16 维系统 (15 eq + 1 arc constraint)
    % branch_follow2 会自行把 [x;mu] 拼成 16 维，所以这里只传 15 个系数
    [x, ~] = branch_follow2('nondim_temp2', branch_len, mu0, mu1, x0(1:15), x0(1:15), sysP);
    
    % 结束后清空全局频率
    FixedOmega = []; 
end