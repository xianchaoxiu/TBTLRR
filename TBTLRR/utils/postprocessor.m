function [W] = postprocessor(Zn)

[~,~,n3] = size(Zn);  % Get the third dimension size
Zn = fft(Zn,[],3);    % Apply FFT along the third dimension

for i = 1:n3
    [U(:,:,i), S, V(:,:,i)] = svd(Zn(:,:,i), 'econ');  % SVD for each slice
    s = diag(S);  % Extract singular values
    r = sum(s > 1e-6);  % Rank approximation

    % Compute matrix M
    U_r = U(:,1:r,i);  % Truncated U
    S_r = diag(s(1:r)); % Truncated singular values
    M = U_r * S_r.^(1/2);

    % Normalize rows of M
    mm = M ./ vecnorm(M, 2, 2); % Normalize each row

    % Compute the weighting matrix W
    rs = mm * mm';
    W(:,:,i) = rs.^2;  % Squaring each element
end
W=ifft(W,[],3);

end
