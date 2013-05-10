function [F,tmpMat]=createReads6Regions_SE(sel,freq,readLength,seq_readLength,values)


min_freq = min(freq);
for i=1:6
    reg = [];
    
    FREQ = [];
    k = 1;
    for j=sel
      curr_freq = round(freq(k)/min_freq*100);
      if seq_readLength(j,(i-1)*readLength+1)~=char(1)
        add_reads_f = seq_readLength(j,(i-1)*readLength+1:i*readLength);
        add_reads_r = seq_readLength(j,(i-1)*readLength+1+readLength*6:i*readLength+readLength*6);
        FREQ = [FREQ;curr_freq;curr_freq];
        reg = [reg;add_reads_f;add_reads_r];        
        
      end
      k = k+1;
    end
%keyboard    
    red = pack_seqs(reg,64);

    [uniqueReads,uniqueReads_inds] = extract_sub_kmers(red, readLength*ones(size(red,1),1),readLength, 1,0);
    clear red
    
    [junk_vals junk_inds tmp_uniqueReads_length] = get_duplicates(uniqueReads_inds(:,1));
    uniqueReads_length = tmp_uniqueReads_length.*(cellfun(@(x) sum(FREQ(uniqueReads_inds(x,2))),junk_inds))';
    clear junk_vals junk_inds uniqueReads_inds

    reads{i}.uniqueReads = uniqueReads;
    reads{i}.uniqueReads_length = uniqueReads_length;
end

for i=1:6
  
    [fracRelevantReads{i},sumRelevantReads{i}] = currReadsFourth(reads{i}.uniqueReads,reads{i}.uniqueReads_length,values{i},0,struct);
    
    F{i} = fracRelevantReads{i}./sum(fracRelevantReads{i});
    
    tmpMat{i} = sparse(size(fracRelevantReads{i},1),6);
    tmpMat{i}(:,i) = -F{i};
end
