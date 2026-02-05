function z = nondim_temp2(y, sysP)
% nondim_temp2
% -------------------------------------------------------------------------
% 兼容两种调用方式：
%   (A) FRF 初值/Newton：y为16维 [15状态; W]，激励力用全局 Fw
%   (B) 弧长连续：y为15维状态，连续参数用全局 mu
%       - FixedOmega 为空：mu 视为频率 W（FRF 连续）
%       - FixedOmega 非空：mu 视为激励力 F（Level-1 / 分岔连续）
%
% 输出始终为 15x1 残差向量（HBM 的 3 组 5 方程：主/副/电路）
% -------------------------------------------------------------------------

    global Fw FixedOmega mu

    len_y = length(y);

    % 解析连续参数（频率W 与 力F）
    if len_y == 16
        % Newton：根据 FixedOmega 判别最后一维是“频率”还是“力”
        state = y(1:15);
        if isempty(FixedOmega)
            % FRF 模式：最后一维是频率
            W = y(16);
            current_Fw = Fw;
        else
            % Level-1 模式：频率固定，最后一维是力
            W = FixedOmega;
            current_Fw = y(16);
        end
    elseif len_y == 15
        % Continuation：状态15维，参数来自全局 mu
        state = y(1:15);
        if isempty(FixedOmega)
            % FRF continuation：mu 作为频率
            W = mu;
            current_Fw = Fw;
        else
            % Level-1 continuation：mu 作为力，频率固定
            W = FixedOmega;
            current_Fw = mu;
        end
    else
        error('nondim_temp2: 维度异常 (%d)。期望 15(continuation) 或 16(Newton-FRF)。', len_y);
    end

    % 提取物理参数
    be1=sysP(1); be2=sysP(2); mu_mass=sysP(3); al1=sysP(4);
    ga1=sysP(5); ze1=sysP(6); lam=sysP(7); sig=sysP(8); kap_e=sysP(9);

    % HBM 系数分组（每组5个：A0, A1c, B1s, A3c, B3s）
    x1 = state(1:5);
    x2 = state(6:10);
    q  = state(11:15);

    % 残差
    R1 = hbm_op(x1, x2, W, 'primary',   be1, ze1, al1, ga1);
    R1(2) = R1(2) - current_Fw; % 外激励作用在主质量 1Ω 余弦项
    R2 = hbm_op(x2, x1, W, 'secondary', be1, ze1, be2, mu_mass, lam, q);
    R3 = hbm_op(q,  x2, W, 'circuit',   sig, kap_e);

    z = [R1; R2; R3];
end

% -------------------------------------------------------------------------
% hbm_op：一阶+三次谐波(1Ω,3Ω) 的代数方程残差
% -------------------------------------------------------------------------
function Res = hbm_op(u, v, W, type, p1, p2, p3, p4, p5, v2)
    Res = zeros(5,1);
    switch type
        case 'primary'
            ga1 = p4; al1 = p3; be1 = p1; ze1 = p2;
            Res(1) = be1*(u(1)-v(1)) + al1*u(1) + ga1*(u(1)^3 + 1.5*u(1)*(u(2)^2+u(3)^2+u(4)^2+u(5)^2));
            Res(2) = -W^2*u(2) + 2*ze1*W*(u(3)-v(3)) + be1*(u(2)-v(2)) + al1*u(2) + ga1*(3*u(1)^2*u(2) + 0.75*u(2)*(u(2)^2+u(3)^2+2*u(4)^2+2*u(5)^2));
            Res(3) = -W^2*u(3) - 2*ze1*W*(u(2)-v(2)) + be1*(u(3)-v(3)) + al1*u(3) + ga1*(3*u(1)^2*u(3) + 0.75*u(3)*(u(2)^2+u(3)^2+2*u(4)^2+2*u(5)^2));
            Res(4) = -9*W^2*u(4) + 6*ze1*W*(u(5)-v(5)) + be1*(u(4)-v(4)) + al1*u(4) + ga1*(3*u(1)^2*u(4) + 0.75*u(4)*(2*u(2)^2+2*u(3)^2+u(4)^2+u(5)^2));
            Res(5) = -9*W^2*u(5) - 6*ze1*W*(u(4)-v(4)) + be1*(u(5)-v(5)) + al1*u(5) + ga1*(3*u(1)^2*u(5) + 0.75*u(5)*(2*u(2)^2+2*u(3)^2+u(4)^2+u(5)^2));

        case 'secondary'
            mu_mass = p4; lam = p5; be1 = p1; be2 = p3; ze1 = p2; q = v2;
            Res(1) = be1*(u(1)-v(1)) + be2*u(1);
            Res(2) = -mu_mass*W^2*u(2) + 2*ze1*W*(u(3)-v(3)) + be1*(u(2)-v(2)) + be2*u(2) + lam*W*q(3);
            Res(3) = -mu_mass*W^2*u(3) - 2*ze1*W*(u(2)-v(2)) + be1*(u(3)-v(3)) + be2*u(3) - lam*W*q(2);
            Res(4) = -9*mu_mass*W^2*u(4) + 6*ze1*W*(u(5)-v(5)) + be1*(u(4)-v(4)) + be2*u(4) + 3*lam*W*q(5);
            Res(5) = -9*mu_mass*W^2*u(5) - 6*ze1*W*(u(4)-v(4)) + be1*(u(5)-v(5)) + be2*u(5) - 3*lam*W*q(4);

        case 'circuit'
            sig = p1; kap_e = p2; xi2 = v;
            Res(1) = kap_e*u(1);
            Res(2) = -W^2*u(2) + sig*W*u(3) + kap_e*u(2) - W*xi2(3);
            Res(3) = -W^2*u(3) - sig*W*u(2) + kap_e*u(3) + W*xi2(2);
            Res(4) = -9*W^2*u(4) + 3*sig*W*u(5) + kap_e*u(4) - 3*W*xi2(5);
            Res(5) = -9*W^2*u(5) - 3*sig*W*u(4) + kap_e*u(5) + 3*W*xi2(4);
    end
end
