
function [ inv_a ] = transformer_t_inverse(A,Phi)

[~,n2,n3]=size(A);

A = ttrans(A, Phi);
inv_a = zeros(n2,n2,n3);
for i = 1 : n3
   inv_a(:,:,i) =  (A(:,:,i)'*A(:,:,i) + eye(n2))\eye(n2);
end
inv_a = ttrans(inv_a, Phi');
end