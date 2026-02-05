clear; clc; close all;
%% HBM vs ODE45 硬验证专用程序
% 功能：
% 1. 验证 HBM 代码计算出的点，在 ODE45 时域仿真中是否真实存在。
% 2. 修正之前的“残差图”错误，计算真实的代数方程残差。
% 3. 生成可直接放入论文的“Method Validation”图表。

%% 1. 参数设置 (需与 main_function 保持一致)
% ---------------------------------------------------------
mu     = 0.1;
al1    = 0.01;
ga1    = 100.0;     % 非线性强度 (注意：你的代码中非线性仅加在主质量上)
ze1    = 0.005;
be1    = 1.0;
be2    = 1.0;
lam    = 0.5;
sig    = 0.1;
kap_e  = 1.0;
sysP = [be1, be2, mu, al1, ga1, ze1, lam, sig, kap_e];

% --- 选择验证点 ---
% 建议选择 S 形曲线下分支的一个稳定点 (误差最小，最容易验证)
Test_Omega = 5.33;  % 验证频率
Test_Force = 2.0;   % 验证激励力

fprintf('===============================================\n');
fprintf('开始执行硬验证 (Validation Process)\n');
fprintf('目标: Omega = %.2f, Force = %.2f\n', Test_Omega, Test_Force);
fprintf('===============================================\n');

%% 2. HBM 求解 (Frequency Domain)
fprintf('[Step 1] 正在计算 HBM 理论解 (1阶+3阶谐波)...\n');

% 构造初值 (猜一个较小的初始值，确保收敛到下分支)
y_guess = zeros(15, 1); 
y_guess(2) = 0.1; % 给主质量基波幅值一点初值

% 定义匿名函数用于 fsolve，绑定参数
% 注意：我们需要欺骗 nondim_temp2，让它以为是在做 FRF 计算
global Fw FixedOmega mu
Fw = Test_Force; 
FixedOmega = Test_Omega; 
mu = Test_Force; % 确保参数传递万无一失

% 调用 fsolve 解代数方程
options = optimoptions('fsolve','Display','off', 'FunctionTolerance', 1e-10);
[y_hbm, fval, exitflag] = fsolve(@(y) nondim_temp2_wrapper(y, sysP), y_guess, options);

% 计算真实的代数残差 (True Residual)
Res_Norm = norm(fval);
Amp_HBM = sqrt(y_hbm(2)^2 + y_hbm(3)^2); % 主质量基波幅值

fprintf('   -> HBM 求解结束. ExitFlag = %d\n', exitflag);
fprintf('   -> 真实代数残差 ||R|| = %.4e (之前看到的 10^-1 是假的)\n', Res_Norm);
fprintf('   -> HBM 预测幅值 A1  = %.6f\n', Amp_HBM);

%% 3. ODE45 仿真 (Time Domain)
fprintf('[Step 2] 正在运行 ODE45 时域积分 (可能需要几秒)...\n');

tspan = [0 1500]; % 跑足够长的时间，确保瞬态消失
y0 = zeros(6, 1); % 初始条件 [x1, x1d, x2, x2d, q, qd] 设为0验证下分支

opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t, y_time] = ode45(@(t,y) sys_ode(t, y, sysP, Test_Omega, Test_Force), tspan, y0, opts);

% 截取最后 20 个周期的数据
T = 2*pi/Test_Omega;
idx_steady = t > (t(end) - 20*T);
t_ss = t(idx_steady);
x1_ss = y_time(idx_steady, 1);

% 计算时域幅值 (峰峰值的一半)
Amp_ODE = (max(x1_ss) - min(x1_ss)) / 2;
fprintf('   -> ODE45 仿真结束.\n');
fprintf('   -> ODE45 稳态幅值   = %.6f\n', Amp_ODE);

%% 4. 误差分析与绘图
Error = abs(Amp_ODE - Amp_HBM) / Amp_HBM * 100;
fprintf('-----------------------------------------------\n');
fprintf('验证结论: 相对误差 = %.2f%%\n', Error);
if Error < 1.0
    fprintf('STATUS: [PERFECT MATCH] 算法准确无误。\n');
else
    fprintf('STATUS: [CHECK NEEDED] 误差稍大，可能是高次谐波影响。\n');
end
fprintf('===============================================\n');

% --- 画图 (可以直接放进论文) ---
figure(200); clf;
set(gcf, 'Color', 'w');

% 子图1: 时域波形对比
subplot(2,1,1);
plot(t_ss, x1_ss, 'b-', 'LineWidth', 1.5); hold on;
% 绘制 HBM 预测的包络线 (幅值)
yline(Amp_HBM, 'r--', 'LineWidth', 2);
yline(-Amp_HBM, 'r--', 'LineWidth', 2);
xlim([t_ss(1), t_ss(1)+5*T]); % 只看前5个周期细节
xlabel('Time \tau'); ylabel('Displacement \xi_1');
title(['Validation: Time History (F=', num2str(Test_Force), ', \Omega=', num2str(Test_Omega), ')']);
legend('ODE45 (Exact)', 'HBM Envelope (Approx)', 'Location', 'best');
grid on; set(gca, 'FontSize', 12);

% 子图2: 细节重构对比 (把 HBM 的谐波还原成波形)
subplot(2,1,2);
% 重构 HBM 波形: x = A0 + A1c*cos + B1s*sin + A3c*cos3 + B3s*sin3
phase = Test_Omega * t_ss;
x1_rec = y_hbm(1) + ...
         y_hbm(2)*cos(phase) + y_hbm(3)*sin(phase) + ...
         y_hbm(4)*cos(3*phase) + y_hbm(5)*sin(3*phase);

plot(t_ss, x1_ss, 'b', 'LineWidth', 2.5); hold on;
plot(t_ss, x1_rec, 'r--', 'LineWidth', 2);
xlim([t_ss(1), t_ss(1)+3*T]); % 放大看细节
xlabel('Time \tau'); ylabel('Displacement \xi_1');
title('Detail Comparison: ODE45 vs HBM Reconstruction');
legend('ODE45', 'HBM (1st+3rd Harmonics)');
grid on; set(gca, 'FontSize', 12);

%% --- 内部函数 (根据你的 nondim_temp2 反推的 ODE) ---
function dydt = sys_ode(t, y, sysP, W, F)
    % 解包参数
    be1=sysP(1); be2=sysP(2); mu_val=sysP(3); al1=sysP(4);
    ga1=sysP(5); ze1=sysP(6); lam=sysP(7); sig=sysP(8); kap_e=sysP(9);
    
    % 状态变量: y = [x1; x1_dot; x2; x2_dot; q; q_dot]
    x1 = y(1); x1d = y(2);
    x2 = y(3); x2d = y(4);
    q  = y(5); qd  = y(6);
    
    dydt = zeros(6,1);
    
    % 1. 主质量方程 (Primary Mass)
    % x1'' + 2*ze1*(x1'-x2') + be1*(x1-x2) + al1*x1 + ga1*x1^3 = F*cos(Wt)
    % 注意：你的 nondim_temp2 里 ga1 只作用于 x1^3 (u(1)^3)，没有耦合非线性
    dydt(1) = x1d;
    dydt(2) = F*cos(W*t) - (2*ze1*(x1d-x2d) + be1*(x1-x2) + al1*x1 + ga1*x1^3);
    
    % 2. 副质量方程 (Secondary Mass)
    % mu*x2'' + 2*ze1*(x2'-x1') + be1*(x2-x1) + be2*x2 + lam*q' = 0
    dydt(3) = x2d;
    dydt(4) = ( - (2*ze1*(x2d-x1d) + be1*(x2-x1) + be2*x2 + lam*qd) ) / mu_val;
    
    % 3. 电路方程 (Circuit)
    % q'' + sig*q' + kap_e*q - x2' = 0
    dydt(5) = qd;
    dydt(6) = x2d - sig*qd - kap_e*q;
end

function R = nondim_temp2_wrapper(y, sysP)
    % 包装器：为了适配 fsolve 输入格式
    % y 是 15 维， nondim_temp2 需要根据全局变量判断
    R = nondim_temp2(y, sysP);
end