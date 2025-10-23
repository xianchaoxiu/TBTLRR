function [Z_star, weights, diag_ratios] = diagonal_ratio_weighted_fusion(Z)
% 基于对角占比的自适应通道加权融合
% 输入：
%   Z: N x N x M 的张量（通道在第三维）
% 输出：
%   Z_star: N x N 加权融合后的矩阵
%   weights: 1 x M 权重
%   diag_ratios: 1 x M 每通道对角占比

[N, ~, M] = size(Z);
weights = zeros(1, M);
diag_ratios = zeros(1, M);

for i = 1:M
    A = Z(:, :, i);
    diag_sum = sum(abs(diag(A)));         % 主对角线元素绝对值和
    total_sum = sum(abs(A(:))) + eps;     % 所有元素绝对值和（加 eps 避免除零）
    diag_ratios(i) = diag_sum / total_sum;
    weights(i) = diag_ratios(i);          % 权重 = 对角占比
end

% 归一化权重
weights = weights / sum(weights);

% 加权融合
Z_star = zeros(N, N);
for i = 1:M
    Z_star = Z_star + weights(i) * Z(:, :, i);
end
end
