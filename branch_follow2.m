function [x, p] = branch_follow2(fname, nsteps, mu0, mu1, x0, x1, sysP)
    %% 优化版弧长连续法
    
    global mu tracking_file_name xc arc
    
    tracking_file_name = fname;
    
    % 初始化
    x0 = x0(:);
    x1 = x1(:);
    
    if numel(x0) == 16
        mu0 = x0(end);
    else
        x0 = [x0; mu0];
    end
    
    if numel(x1) == 16
        mu1 = x1(end);
        xc = x1;
    else
        xc = [x1; mu1];
    end
    
    x = [x0, xc];
    arc = norm(x0 - xc);
    arc_init = arc;
    
    % 自适应参数
    arc_min = arc_init * 0.01;
    arc_max = arc_init * 10;
    res_tol = 1e-5;
    
    k = 1;
    p = 'y';
    success_count = 0;
    fail_count = 0;
    
    fprintf('Branch Follow 开始: 步数=%d, 初始弧长=%.4e\n', nsteps, arc);
    
    while k < nsteps
        % 预测
        if k == 1
            xg = 2*xc - x0;
        else
            % 使用更高阶预测
            if size(x, 2) >= 3
                % 二阶外推
                xg = 3*x(:, end) - 3*x(:, end-1) + x(:, end-2);
            else
                xg = 2*xc - x0;
            end
        end
        
        % 校正
        [xx, ok, Rn] = newton('branch_aux2', xg, sysP);
        
        if ok && Rn < res_tol
            % 成功
            k = k + 1;
            success_count = success_count + 1;
            fail_count = 0;
            
            x0 = xc;
            xc = xx;
            x = [x, xx];
            
            % 自适应调整弧长
            if success_count >= 5 && arc < arc_max
                arc = min(arc * 1.5, arc_max);
                success_count = 0;
            end
            
            % 进度显示
            if mod(k, 50) == 0
                fprintf('  步骤 %d/%d, 参数=%.4f, 残差=%.2e, 弧长=%.2e\n', ...
                        k, nsteps, xc(end), Rn, arc);
            end
            
            % 物理范围检查
            if xc(end) > 300 || xc(end) < 0
                fprintf('参数超出范围，停止。\n');
                break;
            end
            
        else
            % 失败：缩小弧长
            fail_count = fail_count + 1;
            arc = arc * 0.5;
            
            if arc < arc_min
                fprintf('弧长过小(%.2e)，无法继续。\n', arc);
                p = 'n';
                break;
            end
            
            if fail_count > 10
                fprintf('连续失败10次，停止。\n');
                p = 'n';
                break;
            end
        end
    end
    
    fprintf('Branch Follow 完成: 共%d步\n', k);
end