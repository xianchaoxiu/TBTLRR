function X_noisy = add_noise_to_tensor(X, NR, type)
% 添加噪音到张量 X，支持 sparse / saltpepper / gaussian / none
% NR 为噪音比例或标准差（视类型而定）

X_noisy = X;

switch lower(type)
    case 'sparse'
        mask = rand(size(X)) < NR;
        noise = rand(size(X)) .* mask;
        X_noisy = X + noise;

    case 'saltpepper'
        mask = rand(size(X));
        X_noisy(mask < NR/2) = 0;     % 椒
        X_noisy(mask > 1 - NR/2) = 1; % 盐

    case 'gaussian'
        sigma = NR;  % 高斯标准差
        noise = sigma * randn(size(X));
        X_noisy = X + noise;

    case 'none'
        % 无噪音
        X_noisy = X;

    otherwise
        error('Unsupported noise type: %s', type);
end

% 强制截断到 [0,1] 范围（图像场景）
X_noisy = min(max(X_noisy, 0), 1);

end
