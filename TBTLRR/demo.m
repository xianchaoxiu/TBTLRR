clc; clear; close all;
addpath(genpath(cd));
%parpool(6)
%% step 1: 读取数据并处理
DATA=1;
load('extendyaleb.mat') 
maxp=max(fea(:));
if maxp>1
    X=fea/255.0;
else
    X=fea;
end
[n1,n2,n3]=size(X);
if DATA
   %% U matrix
   O = tenmat(X, 3); % unfolding
   O = O.data;
   [U0, D0, V0] = svd(O, 'econ');
   Phi  = U0';
end
%% 添加噪音
NR = 0.35; % sparse noise level
type='sparse';
X = add_noise_to_tensor(X, NR, type);
%% 设置RPCA
opts.denoising_flag = 1; 
% (1 denotes we use R-TPCA; 0 deonotes we do not use)
if opts.denoising_flag
    [n1, n2, n3] = size(X);
    opts.lambda = 5/sqrt(max(n1,n2)*n3);
    opts.mu = 1e-4;
    opts.tol = 1e-4;
    opts.rho = 1.2;
    opts.max_iter = 800;
    opts.DEBUG = 1; 
    opts.Phi = Phi; 
end    
[LL, V, U, RR] = dictionary_learning2(X, opts,Phi);
%% step 2: 设置参数
alpha = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1,10,100];
gamma = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1,10,100];
alpha=1;
gamma=1;
t=0.1;
for a = 1:length(alpha) 
    for b=1:length(gamma)
        tStart = tic;
         [  Z,L,E,J_rank,err ,N,objval] = TBTLRR(X,LL,RR,@soft,alpha(a),gamma(b),Phi);
         Z = ttprod(V, Z,Phi);
         L = ttprod(L, conj_tran(U, Phi),Phi);
 
        elapsedTime = toc(tStart);
        [mean_nmi,BDI,results] = ncut_clustering(Z, gnd,t);
        std1=std(results);
        fprintf('gamma=%.4f  lambda=%.4f  acc=%.4f±%.4f  nmi=%.4f±%.4f  pur=%.4f±%.4f\n', ...
                gamma(b), alpha(a), mean_nmi(1), std1(1), mean_nmi(2), std1(2), mean_nmi(3), std1(3));  
    end
end


