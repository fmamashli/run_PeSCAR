function [unweighted_cond1,weighted_cond1]=compute_permutation_pvalue(sig_matrix_cond1_Allfreq,sig_matrix_cond1_Allfreq_p,nPerm_s)




otc=sum(sum(sig_matrix_cond1_Allfreq~=0));
wotc=sum(sum(sig_matrix_cond1_Allfreq));

count=1;
wcount=1;
for iperm=1:nPerm_s
    
    if sum(sum(squeeze(sig_matrix_cond1_Allfreq_p(iperm,:,:))~=0))>otc
        count=count+1;
    end
    
    if sum(sum(squeeze(sig_matrix_cond1_Allfreq_p(iperm,:,:))))>wotc
        wcount=wcount+1;
    end
    
end

unweighted_cond1=count/(nPerm_s+1);
weighted_cond1=wcount/(nPerm_s+1);