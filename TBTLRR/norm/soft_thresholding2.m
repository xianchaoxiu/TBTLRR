function y = soft_thresholding2(sigma, lambda)
    y = zeros(size(sigma));
    c = nthroot(54, 3) / 4 * lambda^(2/3);  % 阈值

    for i = 1:length(sigma)
        s = sigma(i);
        abs_s = abs(s);

        if abs_s > c
            % 计算 threshold
            phi = (lambda / 8) * (abs_s / 3)^(-2/3);
            phi = min(max(phi, -1), 1);  % 防止 acos 出复数

            theta = acos(phi);
            y(i) = (2/3) * s * (1 + cos(2*pi/3 - (2/3)*theta));
        else
            y(i) = 0;
        end
    end
end
