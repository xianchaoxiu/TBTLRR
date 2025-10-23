function [mean_nmi,BDI,results] = ncut_clustering(Z, labels,t)
warning off;
K = length(unique(labels));
%[Z]=refinecoefficientTensor(Z,ratio);

[n1,n2,n3]=size(Z);
Z=abs(Z);
Z=(Z+tran(Z))./2.0;
%{
XX=zeros(n1,n2);
for i=1:n3
    XX=XX+0.5*((Z(:,:,i)+Z(:,:,i)'));
end
%}
[XX, weights, diag_ratios] = diagonal_ratio_weighted_fusion(Z);

BDI=compute_BDI(XX,labels);

[~,S,U] = svd(XX,'econ');
S = diag(S);
r = sum(S>1e-5*S(1));
U = U(:,1:r);S = S(1:r);
U = U*diag(sqrt(S));
U = normr(U);
L = (U*U').^t;
%% spectral clustering
D = diag(1./sqrt(sum(L,2)));
L = D*L*D;
if any(any(imag(L) ~= 0))
    L1=real(L);
 else
    L1=L;
end

d=50;
results=zeros(50,3);
for p=1:d
    grps = SpectralClustering(L1,K);
    [result] = ClusteringMeasure(labels,grps);
    results(p,:)=result;
end
mean_nmi=mean(results);
end

