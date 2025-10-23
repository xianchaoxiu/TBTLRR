function [X,tnn,trank] = prox_tnn_fun(Y,rho,fun,Phi)
% 自动判断实数/复数变换类型，选择合适的 TNN 近端映射方式

if isreal(Phi)
    % 实数变换（如 DCT/DWT）：使用简洁版本，无需共轭对称
    [X, tnn, trank] = prox_tnn_real(Y, rho, fun, Phi);
else
    % 复数变换（如 DFT）：使用完整版本，保留共轭对称结构
    [X, tnn, trank] = prox_tnn_complex(Y, rho, fun, Phi);
end
end




function [X, tnn, trank] = prox_tnn_real(Y, rho, fun, Phi)

[n1, n2, n3] = size(Y);
X = zeros(n1, n2, n3);
Y = ttrans(Y, Phi);
tnn = 0;
trank = 0;

for i = 1 : n3
    [U, S, V] = svd(Y(:,:,i), 'econ');
    s = diag(S);
    s_thresh = fun(s, rho);
    r = sum(s_thresh > 0);
    if r >= 1
        X(:,:,i) = U(:,1:r) * diag(s_thresh(1:r)) * V(:,1:r)';
        tnn = tnn + (sum(s_thresh.^(1/2)))^2;
        trank = max(trank, r);
    end
end

tnn = tnn / n3;
X = ttrans(X, Phi');
end



function [X,tnn,trank] = prox_tnn_complex(Y,rho,fun,Phi)
[n1,n2,n3] = size(Y);
X = zeros(n1,n2,n3);
Y = ttrans(Y, Phi);
tnn = 0;
trank = 0;
        
% first frontal slice
[U,S,V] = svd(Y(:,:,1),'econ');
S = diag(S);
S0=fun(S,rho);
     r = length(find(S0>0));
    %S = S(1:r)-rho;
if r>=1
    X(:,:,1) = U(:,1:r)*diag(S0(1:r))*V(:,1:r)';
    tnn = tnn+(sum(S0.^(1/2)))^2;
    trank = max(trank,r);
end
% i=2,...,halfn3
halfn3 = round(n3/2);
for i = 2 : halfn3
    [U,S,V] = svd(Y(:,:,i),'econ');
    S = diag(S);
    S0=fun(S,rho);
     r = length(find(S0>0));
    %         S = S(1:r)-rho;
    if r>=1
        X(:,:,i) = U(:,1:r)*diag(S0(1:r))*V(:,1:r)';
        tnn = tnn+(sum(S0.^(1/2)))^2;
        trank = max(trank,r);
    end
    X(:,:,n3+2-i) = conj(X(:,:,i));
end

% if n3 is even
if mod(n3,2) == 0
    i = halfn3+1;
    [U,S,V] = svd(Y(:,:,i),'econ');
    S = diag(S);
      S0=fun(S,rho);
     r = length(find(S0>0));
    %  S = S(1:r)-rho;
    if r>=1
       
        X(:,:,i) = U(:,1:r)*diag(S0(1:r))*V(:,1:r)';
        tnn = tnn+(sum(S0.^(1/2)))^2;
        trank = max(trank,r);
    end
end
tnn = tnn/n3;
X = ttrans(X, Phi');
end
