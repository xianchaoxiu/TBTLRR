function [Phi, Phi_inv] = generate_dwt_level1_basis(n, wname)
% 一层正交小波变换矩阵 Phi 和近似逆 Phi_inv
% 输入:
%   n —— 信号长度（必须是 2 的倍数）
%   wname —— 正交小波名称，如 'haar', 'db1', 'sym4'
%
% 输出:
%   Phi —— 正交小波变换矩阵（一层 DWT）
%   Phi_inv —— 逆变换矩阵，若正交则 Phi_inv ≈ Phi'

    if mod(n, 2) ~= 0
        error('n must be even for level-1 DWT.');
    end

    Phi = zeros(n);
    Phi_inv = zeros(n);

    for i = 1:n
        e = zeros(n, 1);
        e(i) = 1;
        [c, l] = wavedec(e, 1, wname);  % 一层分解
        Phi(:, i) = c(:);               % 保存变换系数
    end

for i = 1:n
    e = zeros(n, 1);
    e(i) = 1;
    [c, l] = wavedec(e, 1, wname);  % 正向获取结构信息
    x_rec = waverec(c, l, wname);   % 使用结构做逆变换
    Phi_inv(i, :) = x_rec(:)';
end


    % 可选：验证是否正交
    ortho_err = norm(Phi_inv * Phi - eye(n));
    fprintf('[DWT Level 1] Φ′Φ ≈ I，正交误差: %.2e\n', ortho_err);
end
