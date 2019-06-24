function stat=clustterstat2D(G1,G2,cfg)
% Cluster based statitics for 2D

% Group1-observation x Time X Freq
% Group2-observation x Time X Freq
% cfg config file

% stats; cluster pvalues

if ~isfield(cfg,'alpha')
cfg.alpha=0.05;
end

if ~isfield(cfg,'numperm')
cfg.numperm=1000;
end

if ~isfield(cfg,'statmethod')
cfg.statmethod='ranksum';
end

if ~isfield(cfg,'conn')
cfg.conn=1;
end


[cluster_P_obs, cluster_N_obs, L_P, L_N,z_obs]=makecluster(G1,G2,cfg);



if isempty(cluster_P_obs)
    display('No Positive Cluster')
else
    [cluster_P_obs, IndexClusP]=sort(cluster_P_obs,'descend');
end
if isempty(cluster_N_obs)
    display('No negative Cluster')
else
    [cluster_N_obs, IndexClusN]=sort(cluster_N_obs);
end


statclust_p=zeros(1,cfg.numperm);
statclust_n=zeros(1,cfg.numperm);

numSizeG1=size(G1,1);
numSizeG2=size(G2,1);

numSizeG=numSizeG1+numSizeG2;
% matlabpool 4
% %POOL=parpool('local',20);
parfor iperm=1:cfg.numperm
    %RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
    G=cat(1,G1,G2);
    
    G=G(randperm(numSizeG),:,:);
    
    G1_temp=G(1:numSizeG1,:,:);
    G2_temp=G(numSizeG1+1:end,:,:);
    
[cluster_P ,cluster_N]=makecluster(G1_temp,G2_temp,cfg);

if isempty(cluster_P)
    statclust_p(iperm)=0;
else
  statclust_p(iperm)=max(cluster_P);  
end


if isempty(cluster_N)
    statclust_n(iperm)=0;
else
    statclust_n(iperm)=min(cluster_N);  
    
end 

end
%delete(POOL)
%matlabpool close

if isempty(cluster_P_obs)
    stat.posclus=[];
else
    for i=1:length(cluster_P_obs)
    stat.posclus(i).clustermass=cluster_P_obs(i);
    valgreater=sum(statclust_p>=cluster_P_obs(i));
    
    
    stat.posclus(i).pvalue=(valgreater + 1)/(1 + cfg.numperm);
    stat.posclus(i).mask=(L_P==IndexClusP(i));
    
    end
end



if isempty(cluster_N_obs)
    stat.negclus=[];
else
    for i=1:length(cluster_N_obs)
    stat.negclus(i).clustermass=cluster_N_obs(i);
    valgreater=sum(statclust_n<=cluster_N_obs(i));
    
    
    stat.negclus(i).pvalue=(valgreater + 1)/(1 + cfg.numperm);
    stat.negclus(i).mask=(L_N==IndexClusN(i));
    
    end
end


stat.L_P=L_P;

stat.L_N=L_N;

stat.cfg=cfg;
stat.z_obs=z_obs;



function [cluster_P, cluster_N, L_P, L_N,zvalue_obs]=makecluster(Group1,Group2,cfg)


if strcmp(cfg.statmethod,'ttest')
df1=size(Group1,1);
df2=size(Group2,1);
df=df1+df2;

zval=abs(tinv(cfg.alpha/2,df));
[~,~,~,statsttest]=ttest2(Group1,Group2,[],[],[],1);
zvalue_obs=squeeze(statsttest.tstat);
clear statsttest

elseif strcmp(cfg.statmethod,'pairedttest')
df1=size(Group1,1);
df2=size(Group2,1);
if df1~=df2
    error('for paired test size need to be same')
end
df=df1-1;

zval=abs(tinv(cfg.alpha/2,df));
[~,~,~,statsttest]=ttest(Group1-Group2,[],[],[],1);
zvalue_obs=squeeze(statsttest.tstat);
clear statsttest


elseif strcmp(cfg.statmethod,'ranksum')

tempTD=reshape(Group1,size(Group1,1),size(Group1,2)*size(Group2,3));
tempASD=reshape(Group2,size(Group2,1),size(Group2,2)*size(Group2,3));



zval=abs(norminv(1-cfg.alpha/2));
zvalue_obs=zeros(1,size(Group1,2)*size(Group2,3));
for i=1:size(Group1,2)*size(Group2,3)
zvalue_obs(i)=myranksum(tempTD(:,i),tempASD(:,i));
if isnan(myranksum(tempTD(:,i),tempASD(:,i)))
  zvalue_obs(i)=0;  
end
end

zvalue_obs=reshape(zvalue_obs,size(Group1,2),size(Group2,3));

end

zvalue_obsBW_P=zvalue_obs >= zval;

zvalue_obsBW_N=zvalue_obs <= -zval;

zvalue_obs(isnan(zvalue_obs))=0;

zvalue_obsBW_P = bwareaopen(zvalue_obsBW_P, cfg.conn);

zvalue_obsBW_N = bwareaopen(zvalue_obsBW_N, cfg.conn);



[L_P, num_P] = bwlabel(zvalue_obsBW_P);


[L_N, num_N] = bwlabel(zvalue_obsBW_N);



if num_P > 0
cluster_P=zeros(1,num_P);

for i_num_P=1:num_P
    temp=(L_P==i_num_P).*zvalue_obs;
    cluster_P(i_num_P)=sum(temp(:));
end

else
    cluster_P=[];
end

if num_N > 0
cluster_N=zeros(1,num_N);

for i_num_N=1:num_N
    temp=(L_N==i_num_N).*zvalue_obs;
    cluster_N(i_num_N)=sum(temp(:));
end

else
    cluster_N=[];
end





