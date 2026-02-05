function [x, converged, R_norm] = newton(fun, x, sysP)
    %% 优化版 Newton 求解器（无并行）
    
    max_iter = 100;
    tol      = 1e-9;
    h_step   = 1e-6;
    min_step = 1e-10;
    
    converged = false;
    n_vars = length(x);
    
    for k = 1:max_iter
        % 计算残差
        R = feval(fun, x, sysP);
        n_eqs = length(R);
        
        % 判断模式
        if n_vars == n_eqs
            active_indices = 1:n_vars;
        elseif n_vars == n_eqs + 1
            active_indices = 1:n_eqs;
        else
            error('Newton: 维度不匹配');
        end
        
        R_norm = norm(R);
        if R_norm < tol
            converged = true;
            return;
        end
        
        % 计算雅可比
        n_active = length(active_indices);
        J = zeros(n_eqs, n_active);
        
        for i = 1:n_active
            idx = active_indices(i);
            x_temp = x;
            h_adaptive = h_step * max(abs(x(idx)), 1);
            x_temp(idx) = x_temp(idx) + h_adaptive;
            R_temp = feval(fun, x_temp, sysP);
            J(:, i) = (R_temp - R) / h_adaptive;
        end
        
        % 求解
        cond_J = cond(J);
        
        if cond_J > 1e14
            % 正则化
            lambda = 1e-3 * trace(J'*J) / n_active;
            dx = -(J' * J + lambda * eye(n_active)) \ (J' * R);
        elseif cond_J > 1e10
            % SVD
            [U, S, V] = svd(J, 'econ');
            s = diag(S);
            tol_svd = max(size(J)) * eps(max(s));
            r = sum(s > tol_svd);
            if r > 0
                dx = -V(:, 1:r) * ((U(:, 1:r)' * R) ./ s(1:r));
            else
                dx = zeros(n_active, 1);
            end
        else
            % QR分解
            dx = -J \ R;
        end
        
        % 线搜索
        alpha = 1.0;
        c1 = 1e-4;
        rho = 0.5;
        
        for ls = 1:20
            x_trial = x;
            x_trial(active_indices) = x_trial(active_indices) + alpha * dx;
            R_trial = feval(fun, x_trial, sysP);
            R_norm_new = norm(R_trial);
            
            if R_norm_new <= R_norm * (1 - c1 * alpha)
                x = x_trial;
                R_norm = R_norm_new;
                break;
            end
            
            alpha = alpha * rho;
            if alpha < min_step
                x(active_indices) = x(active_indices) + min_step * dx;
                break;
            end
        end
        
        % 收敛检查
        if R_norm < tol || norm(dx) < tol * (1 + norm(x(active_indices)))
            converged = true;
            return;
        end
    end
end