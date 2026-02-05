function z = nondim_temp2(y, sysP)
    %% 优化版：向量化计算HBM残差
    
    global Fw FixedOmega mu
    
    len_y = length(y);
    
    % 解析参数
    if len_y == 16
        state = y(1:15);
        if isempty(FixedOmega)
            W = y(16);
            current_Fw = Fw;
        else
            W = FixedOmega;
            current_Fw = y(16);
        end
    elseif len_y == 15
        state = y(1:15);
        if isempty(FixedOmega)
            W = mu;
            current_Fw = Fw;
        else
            W = FixedOmega;
            current_Fw = mu;
        end
    else
        error('维度异常');
    end
    
    % 提取参数
    be1 = sysP(1);
    be2 = sysP(2);
    mu_mass = sysP(3);
    al1 = sysP(4);
    ga1 = sysP(5);
    ze1 = sysP(6);
    lam = sysP(7);
    sig = sysP(8);
    kap_e = sysP(9);
    
    % 分组系数
    x1 = state(1:5);
    x2 = state(6:10);
    q  = state(11:15);
    
    % 预计算常用项
    W2 = W * W;
    W_ze1_2 = 2 * ze1 * W;
    
    % 优化的HBM算子 - 修正参数传递
    R1 = hbm_op_optimized(x1, x2, W, W2, W_ze1_2, be1, ze1, al1, ga1);
    R1(2) = R1(2) - current_Fw;
    
    R2 = hbm_op_secondary(x2, x1, q, W, W2, W_ze1_2, be1, be2, mu_mass, lam, ze1);
    R3 = hbm_op_circuit(q, x2, W, W2, sig, kap_e);
    
    z = [R1; R2; R3];
end

function Res = hbm_op_optimized(u, v, W, W2, W_ze1_2, be1, ze1, al1, ga1)
    % 向量化计算主质量残差
    Res = zeros(5, 1);
    
    % 预计算
    u_sq_sum = u(2)^2 + u(3)^2 + u(4)^2 + u(5)^2;
    u1_sq = u(1)^2;
    
    % DC项
    Res(1) = be1*(u(1) - v(1)) + al1*u(1) + ga1*(u1_sq*u(1) + 1.5*u(1)*u_sq_sum);
    
    % 基波余弦
    term1 = 3*u1_sq*u(2) + 0.75*u(2)*(u(2)^2 + u(3)^2 + 2*u(4)^2 + 2*u(5)^2);
    Res(2) = -W2*u(2) + W_ze1_2*(u(3) - v(3)) + be1*(u(2) - v(2)) + ...
             al1*u(2) + ga1*term1;
    
    % 基波正弦
    term2 = 3*u1_sq*u(3) + 0.75*u(3)*(u(2)^2 + u(3)^2 + 2*u(4)^2 + 2*u(5)^2);
    Res(3) = -W2*u(3) - W_ze1_2*(u(2) - v(2)) + be1*(u(3) - v(3)) + ...
             al1*u(3) + ga1*term2;
    
    % 三次谐波余弦
    term3 = 3*u1_sq*u(4) + 0.75*u(4)*(2*u(2)^2 + 2*u(3)^2 + u(4)^2 + u(5)^2);
    Res(4) = -9*W2*u(4) + 6*ze1*W*(u(5) - v(5)) + be1*(u(4) - v(4)) + ...
             al1*u(4) + ga1*term3;
    
    % 三次谐波正弦
    term4 = 3*u1_sq*u(5) + 0.75*u(5)*(2*u(2)^2 + 2*u(3)^2 + u(4)^2 + u(5)^2);
    Res(5) = -9*W2*u(5) - 6*ze1*W*(u(4) - v(4)) + be1*(u(5) - v(5)) + ...
             al1*u(5) + ga1*term4;
end

function Res = hbm_op_secondary(u, v, q, W, W2, W_ze1_2, be1, be2, mu_mass, lam, ze1)
    % 副质量残差 - 参数顺序修正
    Res = zeros(5, 1);
    
    % DC
    Res(1) = be1*(u(1) - v(1)) + be2*u(1);
    
    % 基波
    mu_W2 = mu_mass * W2;
    lam_W = lam * W;
    
    Res(2) = -mu_W2*u(2) + W_ze1_2*(u(3) - v(3)) + be1*(u(2) - v(2)) + ...
             be2*u(2) + lam_W*q(3);
    Res(3) = -mu_W2*u(3) - W_ze1_2*(u(2) - v(2)) + be1*(u(3) - v(3)) + ...
             be2*u(3) - lam_W*q(2);
    
    % 三次谐波
    mu_W2_9 = 9 * mu_W2;
    lam_W_3 = 3 * lam_W;
    W_ze1_6 = 6 * ze1 * W;
    
    Res(4) = -mu_W2_9*u(4) + W_ze1_6*(u(5) - v(5)) + be1*(u(4) - v(4)) + ...
             be2*u(4) + lam_W_3*q(5);
    Res(5) = -mu_W2_9*u(5) - W_ze1_6*(u(4) - v(4)) + be1*(u(5) - v(5)) + ...
             be2*u(5) - lam_W_3*q(4);
end

function Res = hbm_op_circuit(u, xi2, W, W2, sig, kap_e)
    % 电路残差
    Res = zeros(5, 1);
    
    sig_W = sig * W;
    
    Res(1) = kap_e * u(1);
    Res(2) = -W2*u(2) + sig_W*u(3) + kap_e*u(2) - W*xi2(3);
    Res(3) = -W2*u(3) - sig_W*u(2) + kap_e*u(3) + W*xi2(2);
    Res(4) = -9*W2*u(4) + 3*sig_W*u(5) + kap_e*u(4) - 3*W*xi2(5);
    Res(5) = -9*W2*u(5) - 3*sig_W*u(4) + kap_e*u(5) + 3*W*xi2(4);
end