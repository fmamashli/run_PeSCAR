function [pos,neg,posclus,negclus]=stats_cluster_eval(STATS,thresh)



neg=0;
pos=0;
posclus=0;
negclus=0;

if ~isempty(STATS)
    if ~isempty(STATS.posclus)
        
        for iclus=1:length(STATS.posclus)
            if STATS.posclus(iclus).pvalue<thresh
              %  posclus=posclus+STATS.posclus(iclus).mask.*STATS.z_obs;  
                posclus=posclus+STATS.posclus(iclus).mask;  
                
                pos=pos+STATS.posclus(iclus).clustermass;
                
            end
        end
    end
    
    
    if ~isempty(STATS.negclus)
        
        for iclus=1:length(STATS.negclus)
            if STATS.negclus(iclus).pvalue<thresh
              %  negclus=negclus+STATS.negclus(iclus).mask.*STATS.z_obs;
                negclus=negclus+STATS.negclus(iclus).mask;
                neg=neg+STATS.negclus(iclus).clustermass;
                
            end
        end
    end
    
end

