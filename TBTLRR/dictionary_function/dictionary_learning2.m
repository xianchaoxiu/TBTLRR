function [LL,V,U,RR] = dictionary_learning2(X,opts,Phi)
if opts.denoising_flag==0
%% directly using raw data as dictionary
    tho=100;% sigma(i)<=tho*sigma(1)
    [ ~,~,U,V,S ] = prox_tran_low_rank_dct(X,tho,Phi);
    LL=ttprod(U,S,Phi);
    RR = ttprod(S,conj_tran(V,Phi),Phi);
else
%% use R-TPCA to denoise data first and then use the recovered data as dictionary
    %% raw R-TPCA algorithm
    [L,~,~,~,~,~] =transformer_trpca_tnn(X,opts.lambda,opts);
    %% approximate L, since sometimes R-TPCA cannot produce a good dictionary
    [U,V,S,~] = tran_tSVDs(L,Phi);
    LL = ttprod(U,S,Phi);
    RR = ttprod(S,conj_tran(V,Phi),Phi);
end
