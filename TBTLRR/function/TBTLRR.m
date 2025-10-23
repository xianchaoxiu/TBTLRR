function [ Z,L,E,J_rank,err ,N,objval] = TBTLRR(X,A,B,fun,alpha,gamma,Phi)

[n1,n2,n3]=size(X);
[~,r,~]=size(A);


%% Z=J=Y2 n2 n2 n3
Z=zeros(r,n2,n3);
J=Z;
Y2=Z;
%% L=S=Y3 n1 n1 n3
L=zeros(n1,r,n3);
S=L;
Y3=L;
%% E=Y1 n1 n2 n3
E = zeros(n1,n2,n3);
Y1=E;
%% N=Y1 n1 n2 n3
N= zeros(n1,n2,n3);
Debug=0;
lambda=alpha;
beta=1e-3;
max_beta = 1e+10;
tol = 1e-7;
rho = 1.5;
iter = 0;
max_iter = 150;
objval = zeros(max_iter, 3); 
%%?Pre compute
AT = conj_tran(A,Phi);
BT = conj_tran(B,Phi);
Ain=transformer_t_inverse(A,Phi);
Bin=transformer_t_inverse2(B,Phi);
Err = zeros(max_iter,3);
while iter < max_iter
    iter = iter+1;
   %% upate Jk
    J_pre = J;
    R1 = Z+Y2/beta;
    [J,J_tnn,J_rank] =prox_tnn_fun(R1,1/beta,fun,Phi);

    %% update Sk
    S_pre = S;
    R2 = L+Y3/beta;
    [S,S_tnn,~] = prox_tnn_fun(R2,1/beta,fun,Phi);
    %% update Zk
    Z_pre=Z;
    Q1=X-ttprod(L,B,Phi)-E-N;
    Z=ttprod(Ain,J+ttprod(AT,Y1,Phi)/beta-Y2/beta+ttprod(AT,Q1,Phi),Phi);
    
    %% update Lk
    L_pre=L;
    Q2=X-ttprod(A,Z,Phi)-E-N;
    L=ttprod(S+ttprod(Y1,BT,Phi)/beta-Y3/beta+ttprod(Q2,BT,Phi),Bin,Phi);
     %% update Ek
    E_pre = E;
    R3=X-ttprod(A,Z,Phi)-ttprod(L,B,Phi)-N+Y1/beta;
    [E] = prox_l0( R3, lambda/beta);
   
   %% update Nk
   N_pre=N;
   R4=Y1+beta*(X-ttprod(A,Z,Phi)-ttprod(L,B,Phi)-E);
   N=R4/(2*gamma+beta);





    %% check convergence
    leq1 = Z-J;
    leq2 = X-ttprod(A,Z,Phi)-ttprod(L,B,Phi)-E-N;
    leq3 = L-S;

    leqm1 = max(abs(leq1(:)));
    leqm2 = max(abs(leq2(:)));
    leqm3 = max(abs(leq3(:)));
    Err(iter,1)=leqm1;
    Err(iter,2)=leqm2;
    Err(iter,3)=leqm3;

    err = max([leqm1,leqm2,leqm3]);
    Y2T=conj_tran(Y2,Phi);
    Y1T=conj_tran(Y1,Phi);
    Y3T=conj_tran(Y3,Phi);

  objval=Err;
    
    if (Debug && (iter==1 || mod(iter,20)==0))
        sparsity=length(find(E~=0));
        objval(iter,1)=J_tnn+S_tnn+lambda*norm(E(:),1)+gamma*norm(N, 'fro')^2;
      %  fprintf('iter = %d, obj = %.3f, err = %.8f, beta=%.2f, gamma = %d, sparsity=%d\n'...
        %    , iter,J_tnn+S_tnn+lambda*norm(E(:),1)+gamma*norm(X, 'fro')^2,err,beta,gamma,sparsity);
       fprintf('iter = %d, obj = %.3f, err = %.8f, Z_tnn=%.4f, S_tnn = %.4f, E=%.4f,N=%.4f\n'...
            , iter,J_tnn+S_tnn+lambda*norm(E(:),1)+gamma*norm(N, 'fro')^2,err,J_tnn,S_tnn,lambda*norm(E(:),1),gamma*norm(N, 'fro')^2);
    end
    if err < tol 
        break;
    end
    
    %% update Lagrange multiplier and  penalty parameter beta
    Y1 = Y1 +beta*(X-ttprod(A,Z,Phi)-ttprod(L,B,Phi)-E-N);
    Y2 = Y2 +beta*(Z-J);
    Y3 = Y3 +beta*(L-S);
    beta = min(beta*rho,max_beta);
end
end
