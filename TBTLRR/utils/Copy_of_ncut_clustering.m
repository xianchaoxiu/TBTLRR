function [] =  Copy_of_ncut_clustering(Z, labels,t) 
warning off;
K = length(unique(labels));
%[Z]=refinecoefficientTensor(Z,ratio);

[n1,n2,n3]=size(Z);
Z=abs(Z);
Z=(Z+tran(Z))./2.0;

XX=zeros(n1,n2);
for i=1:n3
    XX=XX+0.5*((Z(:,:,i)+Z(:,:,i)'));
end

%[XX, weights, diag_ratios] = diagonal_ratio_weighted_fusion(Z);

BDI=compute_BDI(XX,labels);
%%
    XX(XX < 0.05) = 0.01;
    XX(XX > 0.05) = 0.02;
    figure;
    imagesc(XX);
    axis off
    colormap parula
    colorbar
    drawnow
end

