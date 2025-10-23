function [ inv_a ] = transformer_t_inverse2(A,Phi)

[n1,~,n3]=size(A);

A = ttrans(A, Phi);
inv_a = zeros(n1,n1,n3);
for i = 1 : n3
    inv_a(:,:,i) =  (A(:,:,i)*A(:,:,i)' + eye(n1))\eye(n1);
end
inv_a = ttrans(inv_a, Phi');
