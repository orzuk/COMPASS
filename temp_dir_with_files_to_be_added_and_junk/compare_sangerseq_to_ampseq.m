function [align_vec_mat,dist,xi,classSeqs]=compare_sangerseq_to_ampseq(Sequence_amp,n,seq,F_flag,base_opt,s1_start,s2_start,ampHeader,tit)

%len_seq=min([length(seq),500]);
len_seq = length(seq);

n = min(n,len_seq-s1_start-2); n


if F_flag==1
    % n=400;
    align_mat=zeros(length(Sequence_amp),n+1);
%     sumalign=zeros(length(Sequence_amp),1);
    align_vec_mat=zeros(length(Sequence_amp),n+1);
    dist=zeros(length(Sequence_amp),1);

    s1_to_match = cell(1,length(Sequence_amp));
    
    for i=1:length(Sequence_amp)
        %keyboard
        s1=Sequence_amp{i}(s1_start:len_seq);
        s2=seq(s2_start:len_seq); % Sanger major allele 
        AlignStruct = localalign(s1, s2,'Alphabet','NT');
        start=AlignStruct.Start;start
        s1 = s1(start(1):(start(1)+n));
        s2 = s2(start(2):(start(2)+n));
        % AlignStruct = localalign(s1, s2,'Alphabet','NT');
        %     sumalign(i) = sum(Alignment(2,:)=='|');
        tmp_base_opt = base_opt(s2_start:len_seq,:);
        [dist(i),align_vec_mat(i,:)] = seqDistanceWithAmbiguity(s1,tmp_base_opt(start(2):(start(2)+n),:));
        s1_to_match{i} = s1;
        % align_vec_mat is 1 when sequence is similar to Sanger major allele or Sanger minor allele
        align_mat(i,:) = (s1==s2); % matrix of 1 when sequence matches Sanger major allele and 0 otherwise
        
        %         sumalign(i) = seqpdist([{s1};{s2}],'Alphabet','NT','Method','p-distance','PairwiseAlignment',false);
    end
    
    
    % extract
    
    % align_mat=align_vec_mat=0 mean equal to major
    % align_vec_mat=1 & align_mat=0 - equal to minor but not to major in the sequence
    % align_mat=1 & align_vec_mat==0 can not be
    
    if ~isempty(find(align_mat==1 & align_vec_mat==0))
      disp('problem');keyboard
    end
    
    %keyboard
    [~,xi] = sort(dist);
    align_vec_mat = align_vec_mat(xi,:); % sort the rows
    align_mat = align_mat(xi,:); 
    align_vec_mat(align_vec_mat>align_mat) = 2; 
    
    % color code: 
    % black - (0) - sequence different from major allele and of minor allele (if exists)
    % gray - (1) - sequence same as Sanger major allele
    % white - (2) - sequence same as minor allele
    
    
    % find amplified that match the sanger aligned regions
    distMatOfComparedSeqs = squareform(seqpdist(s1_to_match,'Alphabet','NT','Method','p-distance','PairwiseAlignment',false));
    [xx,yy] = find(triu(distMatOfComparedSeqs,1)==0);
    del = find(xx>=yy);
    xx(del) = [];
    yy(del) = [];
     
    classSeqs = unionfind([xx,yy],size(distMatOfComparedSeqs,1));
    classSeqs = classSeqs(xi);
    
    
elseif F_flag==2
disp('add')
return
    align_mat=zeros(length(Sequence_amp),n+1);
    align_vec_mat=zeros(length(Sequence_amp),n+1);
    dist=zeros(length(Sequence_amp),1);
    
    for i=1:length(Sequence_amp)
        %s1=seqrcomplement(Sequence_amp{i}(end-len_seq:end-s1_start));
        s1=seqrcomplement(Sequence_amp{i}(max(1,length(Sequence_amp{i})-len_seq):end-s1_start));
        s2=seq(s2_start:len_seq); % Sanger major allele 
        AlignStruct = localalign(s1, s2,'Alphabet','NT');
        start=AlignStruct.Start;
        s1 = s1(start(1):(start(1)+n));
        s2 = s2(start(2):(start(2)+n));
        %AlignStruct = localalign(s1, s2,'Alphabet','NT');
        tmp_base_opt = base_opt(s2_start:len_seq,:);
        [dist(i),align_vec_mat(i,:)] = seqDistanceWithAmbiguity(s1,tmp_base_opt(start(2):(start(2)+n),:));
        align_mat(i,:) = s1==s2;
%         sumalign(i) = seqpdist([{s1};{s2}],'Alphabet','NT','Method','p-distance','PairwiseAlignment',false);
    end
    
    [~,xi] = sort(dist);
    align_vec_mat = align_vec_mat(xi,:);
    align_mat = align_mat(xi,:);
    align_vec_mat(align_vec_mat>align_mat) = 2;
    
end



