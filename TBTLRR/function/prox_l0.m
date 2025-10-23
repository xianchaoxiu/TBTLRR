function H = prox_l0(Y, lambda)
    % 迭代半阈值算法 (Iterative Half Thresholding Algorithm)
    % 用于三阶张量的 l_{1/2} 正则化
    % 输入:
    %   Y      - 三阶张量 (p × q × n)
    %   lambda - 正则化参数 (λ > 0)

    % 输出:
    %   H      - 半阈值化后的张量 H_lambda(Y)
    
    % 计算阈值
    threshold = (3 * sqrt(54) / 4) * (lambda)^(2/3);
    
    % 遍历张量 Y 的每个元素
    H = Y; % 先复制 Y
    for i = 1:size(Y,1)
        for j = 1:size(Y,2)
            for k = 1:size(Y,3)
                y_val = Y(i,j,k);
                
                % 判断是否满足迭代条件
                if abs(y_val) > threshold
                    % 计算 Psi_lambda
                    Psi_lambda = acos((lambda / 8) * (abs(y_val) / 3)^(-3/2));
                    
                    % 计算 f_lambda_1_2
                    H(i,j,k) = (2/3) * y_val * (1 + cos(2*pi/3 - (2/3) * Psi_lambda));
                else
                    H(i,j,k) = 0; % 低于阈值的直接置零
                end
            end
        end
    end
end
