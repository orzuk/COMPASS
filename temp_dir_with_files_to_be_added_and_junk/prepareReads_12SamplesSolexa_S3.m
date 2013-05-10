

%%% SHAY samples
if 1==1
  % prepare the reads in files - 100
  clear
  load('~/CS/BAC/12Samples/Solexa/illumina_reads_S_samples')
  sampleName = {'S1','S2','S3','S4'};
  %       (Sample order S1,S2,S3,S4)

  readLength = 100;
  clear red uniqueReads uniqueReads_length
  for i=3
    clear red uniqueReads uniqueReads_length
    red = cellseq100_sample{i};
    
    found = [];
    for j=1:size(red,1)
      if ~isempty(find(red(j,:)~='A' & red(j,:)~='C'& red(j,:)~='G'& red(j,:)~='T'))
        found = [found,j];
       
      end
      if mod(j,10^5)==1
        disp(j)
      end
    end
    red(found,:) = [];
    
    red_packed = pack_seqs(red,64);
    clear red
    
    n = size(red_packed,1);
    groupSize = 10^6;
    part = 1:groupSize:n;
    part(end) = n+1;
    clear tmp_* junk*
    tmp_uniqueReads = [];
    tmp_uniqueReads_length = [];
    for j=1:length(part)-1
      tmp_red = red_packed(part(j):part(j+1)-1,:);
      [junk_uniqueReads,junk_uniqueReads_inds] = extract_sub_kmers(tmp_red, readLength*ones(size(tmp_red,1),1),readLength,1,0);
      [junk_junk_vals junk_junk_inds junk_uniqueReads_length] = get_duplicates(junk_uniqueReads_inds(:,1));
      tmp_uniqueReads = [tmp_uniqueReads;junk_uniqueReads];
      tmp_uniqueReads_length = [tmp_uniqueReads_length,junk_uniqueReads_length];
    end
    
    % unique of unique
    [uniqueReads,junk_uniqueReads_inds] = extract_sub_kmers(tmp_uniqueReads, readLength*ones(size(tmp_uniqueReads,1),1),readLength,1,0);
    [junk_junk_vals junk_junk_inds uniqueReads_length] = get_duplicates(junk_uniqueReads_inds(:,1));
    
    clear red red_packed tmp* junk*
    
    save(['~/CS/BAC/12Samples/Solexa/data/readsSolexa_100_sample_',sampleName{i}],'uniqueReads','uniqueReads_length')
  end
end

if 1==1
  % prepare the reads in files - 50
  clear
  load('~/CS/BAC/12Samples/Solexa/illumina_reads_S_samples')
  sampleName = {'S1','S2','S3','S4'};
%       (Sample order S1,S2,S3,S4)

  readLength = 50;
  clear red uniqueReads uniqueReads_length
  for i=[3]
    clear red uniqueReads uniqueReads_length
    red = cellseq50_sample{i};
    
    found = [];
    for j=1:size(red,1)
      if ~isempty(find(red(j,:)~='A' & red(j,:)~='C'& red(j,:)~='G'& red(j,:)~='T'))
        found = [found,j];
       
      end
      if mod(j,10^5)==1
        disp(j)
      end
    end
    red(found,:) = [];
    
    red_packed = pack_seqs(red,64);
    clear red
    
    n = size(red_packed,1);
    groupSize = 10^6;
    part = 1:groupSize:n;
    part(end) = n+1;
    clear tmp_* junk*
    tmp_uniqueReads = [];
    tmp_uniqueReads_length = [];
    for j=1:length(part)-1
      tmp_red = red_packed(part(j):part(j+1)-1,:);
      [junk_uniqueReads,junk_uniqueReads_inds] = extract_sub_kmers(tmp_red, readLength*ones(size(tmp_red,1),1),readLength,1,0);
      [junk_junk_vals junk_junk_inds junk_uniqueReads_length] = get_duplicates(junk_uniqueReads_inds(:,1));
      tmp_uniqueReads = [tmp_uniqueReads;junk_uniqueReads];
      tmp_uniqueReads_length = [tmp_uniqueReads_length,junk_uniqueReads_length];
    end
    
    % unique of unique
    [uniqueReads,junk_uniqueReads_inds] = extract_sub_kmers(tmp_uniqueReads, readLength*ones(size(tmp_uniqueReads,1),1),readLength,1,0);
    [junk_junk_vals junk_junk_inds uniqueReads_length] = get_duplicates(junk_uniqueReads_inds(:,1));
    
    clear red red_packed tmp* junk*
    
    save(['~/CS/BAC/12Samples/Solexa/data/readsSolexa_50_sample_',sampleName{i}],'uniqueReads','uniqueReads_length')
  end
end

% end SHAY samples





    


% test 
%%%%%%%%%%%%%%%%%%%%%%%%
if 1==2

    clear
    load('~/CS/BAC/12Samples/Solexa/illumina_reads_S_samples','cellseq50_sample')
    sampleName = {'S1','S2','S3','S4'};
    %       (Sample order S1,S2,S3,S4)

    readLength = 50;
    clear red uniqueReads uniqueReads_length
  
    i=3;
    clear red uniqueReads uniqueReads_length
    red = cellseq50_sample{i};
    red_packed = pack_seqs(red(1:2*10^6,:),64);
    
  
    n = 5*10^5;
    groupSize = 5*10^4;
    part = 1:groupSize:n;
    part(end) = n+1;
    clear tmp_* junk*
    tmp_uniqueReads = [];
    tmp_uniqueReads_length = [];
    for j=1:length(part)-1
      tmp_red = red_packed(part(j):part(j+1)-1,:);
      
      [junk_uniqueReads,junk_uniqueReads_inds] = extract_sub_kmers(tmp_red, readLength*ones(size(tmp_red,1),1),readLength,1,0);
      
      [junk_junk_vals junk_junk_inds junk_uniqueReads_length] = get_duplicates(junk_uniqueReads_inds(:,1));
      
      tmp_uniqueReads = [tmp_uniqueReads;junk_uniqueReads];
      
      tmp_uniqueReads_length = [tmp_uniqueReads_length,junk_uniqueReads_length];
      
      
    end
    
    % unique of unique
    [junk_uniqueReads_all,junk_uniqueReads_inds] = extract_sub_kmers(tmp_uniqueReads, readLength*ones(size(tmp_uniqueReads,1),1),readLength,1,0);
    [junk_junk_vals junk_junk_inds junk_uniqueReads_length] = get_duplicates(junk_uniqueReads_inds(:,1));
    
    %[j,i1,i2] = intersect(junk_uniqueReads_all,full_uniqueReads(20869,:),'rows');
    
    L = zeros(length(junk_junk_inds),1);
    for j=1:length(junk_junk_inds)
      L(j) = sum(tmp_uniqueReads_length(junk_uniqueReads_inds(junk_junk_inds{j},2)));
    end
    
    
    [full_uniqueReads,full_uniqueReads_inds] = extract_sub_kmers(red_packed(1:n,:), readLength*ones(n,1),readLength,1,0);
    [full_junk_vals full_junk_inds full_uniqueReads_length] = get_duplicates(full_uniqueReads_inds(:,1));
    
    size(setdiff(junk_uniqueReads_all,full_uniqueReads,'rows'))
    size(setdiff(full_uniqueReads,junk_uniqueReads_all,'rows'))
    
    setdiff(L,full_uniqueReads_length)
    setdiff(full_uniqueReads_length,L)
    
end

