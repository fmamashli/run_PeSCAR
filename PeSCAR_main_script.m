clear all
close all


load('example_simulated_data.mat')

time=cond{1,1}.time;
freq=cond{1,1}.freq;
% cond1 data: 8 subjects * N subroi1 * N subroi2 * freq * time
% cond2 data: 8 subjects * N subroi1 * N subroi2 * freq * time

cond1cond2=[cond{1}.data_subj;cond{2}.data_subj];

roi1 = cond{1,1}.roi1;
roi2 = cond{1,1}.roi2;

sig_matrix_cond1_Allfreq_p=zeros(nPerm_s,length(roi1),length(roi2));
sig_matrix_cond2_Allfreq_p=zeros(nPerm_s,length(roi1),length(roi2));

L=size(cond{1}.data_subj,1);

nperm=1000;
statsmethod='pairedttest';
ALPHA=0.05;


nPerm_s=250;
Thresh=0.05;

for iPerm_s=1:nPerm_s
    
    ind_perm=randperm(size(cond1cond2,1));
    
    cond1=cond1cond2(ind_perm(1:L),:,:,:,:);
    cond2=cond1cond2(ind_perm(L+1:end),:,:,:,:);
    
    
    for iLabel1=1:length(roi1)
        
        
        for iLabel2=1:length(roi2)
            
            iPerm_s
            
            
            tic
            
            cond1D=squeeze(cond1(:,iLabel1,iLabel2,:,:));
            cond2D=squeeze(cond2(:,iLabel1,iLabel2,:,:));
            
            
            STATS=do_stats2D(cond1D,cond2D,nperm,statsmethod,ALPHA);
            
            
            mask1=0;mask2=0;
            clussum_pos=0;clussum_neg=0;
            
            
            if ~isempty(STATS)
                [clussum_pos,clussum_neg,mask1,mask2]=stats_cluster_eval(STATS,Thresh);
            end
            
            
            sig_matrix_cond1_Allfreq_p(iPerm_s,iLabel1,iLabel2)=clussum_pos;
            
            sig_matrix_cond2_Allfreq_p(iPerm_s,iLabel1,iLabel2)=clussum_neg;
            
            
            toc
            
            clear STATS
            
            
        end
        
    end
end


%% PeSCAR p-value computation
% load example permuted data

load('permutation_0_0_0_0_nr_1_snr_0.051_templ_tempr_3sub_norand_15to20f_8subj_stg9parts_0_0_0_0_nr_1_snr_0.051_templ_tempr_3sub_norand_15to20f_8subj_stg9parts_nPerm_s500.mat')

[unweighted_cond1,weighted_cond1]=compute_permutation_pvalue(sig_matrix_cond1_Allfreq,sig_matrix_cond1_Allfreq_p,nPerm_s);
[unweighted_cond2,weighted_cond2]=compute_permutation_pvalue(abs(sig_matrix_cond2_Allfreq),abs(sig_matrix_cond2_Allfreq_p),nPerm_s);

% to adjust for the two comparison, multiply each p-value by 2, Bonferroni
% correction

%% Get the Time-frequency map


tag='0_0_0_0_nr_1_snr_0.05_templ_tempr_3sub_norand_15to20f_8subj_stg9parts_0_0_0_0_nr_1_snr_0.05_templ_tempr_3sub_norand_15to20f_8subj_stg9parts';

posclus=0;
negclus=0;

pos=zeros(length(roi1),length(roi2));
neg=zeros(length(roi1),length(roi2));


for iLabel1=1:length(roi1)
    
    for iLabel2=1:length(roi2)
        
        [iLabel1 iLabel2]
        
        
        load(['./coherence/stats_coh_' roi2{iLabel1} '_' roi1{iLabel2}   ...
            '_' tag '.mat']);
        
        if ~isempty(STATS)
            
            if ~isempty(STATS.posclus)
                
                
                for iclus=1:length(STATS.posclus)
                    if STATS.posclus(iclus).pvalue<0.05
                        posclus=posclus+STATS.posclus(iclus).mask;
                        pos(iLabel1,iLabel2)=STATS.posclus(iclus).clustermass;
                        
                    end
                end
            end
            
            
            if ~isempty(STATS.negclus)
                
                for iclus=1:length(STATS.negclus)
                    if STATS.negclus(iclus).pvalue<0.05
                        negclus=negclus+STATS.negclus(iclus).mask;
                        neg(iLabel1,iLabel2)=STATS.negclus(iclus).clustermass;
                        
                    end
                end
            end
            
        end
        
    end
    
end


figure;
subplot(2,2,1)

imagesc(time*1000,freq,posclus);axis xy;colorbar
colormap('jet')

subplot(2,2,2)
imagesc(pos);axis xy;colorbar;title('cond1>cond2-positive clusters')
colormap('jet')

subplot(2,2,3)
imagesc(time*1000,freq,negclus);axis xy;colorbar;
colormap('jet')


subplot(2,2,4)
imagesc(neg);axis xy;colorbar;title('cond2>cond1-negative clusters')
colormap('jet')





