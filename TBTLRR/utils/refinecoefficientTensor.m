function [Z] = refinecoefficientTensor(C, ratio)
    % 获取张量的尺寸
    [m, n, p] = size(C);
    % 初始化输出张量
    Z = zeros(m, n, p);
    
    % 遍历张量的第三维度
    for k = 1:p
        % 提取第三维度切片
        C_slice = C(:,:,k);
        
        % 对每个二维矩阵进行操作
        for i = 1:m
            for j = 1:n
                % 获取当前系数的绝对值
                Crow = abs(C_slice(i,j,:));  % 获取当前系数的绝对值
                sum_Crow = sum(Crow);  % 当前行的总和
                
                % 对系数进行排序并找到排序后的系数
                [sorted, ~] = sort(Crow, "descend");  % 只关心排序后的系数，不需要索引
                
                % 选择系数的逻辑
                t = 0;
                l = 1;
                if ratio < 1
                    % 选择比例小于1时逐步加和直到满足条件
                    while t/sum_Crow < ratio && l <= numel(sorted)
                        t = t + sorted(l);
                        l = l + 1;
                    end
                else
                    % 如果比例大于等于1,直接选择前`ratio`个
                    l = min(round(ratio), numel(sorted));  % 确保l不会超过sorted的长度
                end
                
                % 更新Z张量的相应位置
                Z(i,j,k) = sorted(1:l);  % 赋值给输出张量
            end
        end
    end
end

