function BDI = compute_BDI(W, labels)
    % 计算矩阵 W 的块对角指数（BDI）
    % W:  权重矩阵 (N x N)
    % labels: 每个数据点的簇标签 (N x 1)
    
    % 获取簇的数量
    clusters = unique(labels);
    c = length(clusters);
    
    % 计算 ||W||_l1 (L1 范数, 即所有元素绝对值之和)
    W_L1 = sum(abs(W(:)));
    
    % 初始化 BDI 计算变量
    intra_cluster_sum = 0;
    inter_cluster_sum = 0;
    
    % 遍历每个簇
    for k = 1:c
        % 获取属于当前簇的索引
        Sk = find(labels == clusters(k));
        S_not_k = find(labels ~= clusters(k));
        
        % 计算块对角部分的和 (同一簇内)
        intra_cluster_sum = intra_cluster_sum + sum(abs(W(Sk, Sk)), 'all');
        
        % 计算跨簇的和 (不同簇间)
        inter_cluster_sum = inter_cluster_sum + sum(abs(W(Sk, S_not_k)), 'all');
    end
    
    % 计算 BDI 指数
    BDI = (intra_cluster_sum - 2 * inter_cluster_sum) / W_L1;
end
