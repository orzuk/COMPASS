function [align_vec_mat,align_mat,overlapData,res,inds,SangerOffset,writeErrors]=compare_sangerseq_to_ampseq2(res,lengthValidSanger,seqSanger,F_flag,base_opt,SangerSeq_start,AmplifiedSeq_start)
%keyboard
%len_seq=min([length(seqSanger),500]);
len_seqSanger = length(seqSanger);

%lengthValidSanger = min(lengthValidSanger,len_seqSanger-SangerSeq_start-2); lengthValidSanger

curr_lengthValidSanger = min([lengthValidSanger,len_seqSanger-AmplifiedSeq_start-2]);

% older version
% if F_flag==1
%   keyboard
%   s1 = res{1}{1}.amp.Sequence(AmplifiedSeq_start:len_seqSanger);
%   s2 = seqSanger(SangerSeq_start:len_seqSanger); % Sanger major allele 
%   AlignStruct = localalign(s1, s2,'Alphabet','NT');
%   AlignStruct.Readme = 'the first is amplified';
%   curr_lengthValidSanger = min([length(s1)-AlignStruct.Start(1),length(s2)-AlignStruct.Start(2),lengthValidSanger-AmplifiedSeq_start]);

% else
%   s1 = seqrcomplement(res{1}{1}.amp.Sequence(max(1,length(res{1}{1}.amp.Sequence)-len_seqSanger):end));
%   s2 = seqSanger(SangerSeq_start:len_seqSanger); % Sanger major allele 
%   AlignStruct = localalign(s1, s2,'Alphabet','NT');
%   AlignStruct.Readme = 'the first is amplified';
%   curr_lengthValidSanger = min([length(s1)-AlignStruct.Start(1),length(s2)-AlignStruct.Start(2),lengthValidSanger-AmplifiedSeq_start+1]);
% end
%keyboard


% n=400;
align_mat=zeros(length(res),curr_lengthValidSanger+1);
align_vec_mat=zeros(length(res),curr_lengthValidSanger+1);
dist=zeros(length(res),1);

s1_to_match.seq = cell(1,length(res));

if F_flag==1
  for i=1:length(res)
    %keyboard
    
    %if res{i}{1}.Full.Numbers==220048,keyboard,end 
      
    % amplified
    %s1 = res{i}{1}.amp.Sequence(1:len_seqSanger);
    s1 =  res{i}{1}.amp.Sequence(AmplifiedSeq_start:len_seqSanger);
    s2 = seqSanger(SangerSeq_start:len_seqSanger); % Sanger major allele 
    
    
    AlignStruct = localalign(s1, s2,'Alphabet','NT');
    AlignStruct.Readme = 'the first is amplified';
    
    
    

    
    
    s1 = s1(AlignStruct.Start(1):(AlignStruct.Start(1)+curr_lengthValidSanger)); % s1 is the amplified to later compared to
    s2 = s2(AlignStruct.Start(2):(AlignStruct.Start(2)+curr_lengthValidSanger));
    
    %if AlignStruct.Start(1)<AlignStruct.Start(2)
    %  disp('problem')
    %  keyboard
    %end
   
   
    tmp_base_opt = base_opt(SangerSeq_start:len_seqSanger,:);
    tmp_base_opt = tmp_base_opt(AlignStruct.Start(2):(AlignStruct.Start(2)+curr_lengthValidSanger),:);
    [dist(i),align_vec_mat(i,:)] = seqDistanceWithAmbiguity(s1,tmp_base_opt);
    
    % align_vec_mat is 1 when sequence is similar to Sanger major allele or Sanger minor allele
    align_mat(i,:) = (s1==s2); % matrix of 1 when sequence matches Sanger major allele and 0 otherwise
    
    align_vec_mat(i,align_vec_mat(i,:)>align_mat(i,:)) = 2; % change minor to major
        
    % matched amplified sequence to later compare
    s1_to_match.seq{i} = s1;
    s1_to_match.SangerOffset(i) = SangerSeq_start+AlignStruct.Start(2)-1;
    
        
    writeErrors =  char('P'*ones(3,length(s1)));
    
    posError_A_major = find(align_vec_mat(i,:)==0 & tmp_base_opt(:,1)'==1);
    posError_C_major = find(align_vec_mat(i,:)==0 & tmp_base_opt(:,2)'==1);
    posError_G_major = find(align_vec_mat(i,:)==0 & tmp_base_opt(:,3)'==1);
    posError_T_major = find(align_vec_mat(i,:)==0 & tmp_base_opt(:,4)'==1);
    
    posError_A_minor = find(align_vec_mat(i,:)==0 & tmp_base_opt(:,1)'==2);
    posError_C_minor = find(align_vec_mat(i,:)==0 & tmp_base_opt(:,2)'==2);
    posError_G_minor = find(align_vec_mat(i,:)==0 & tmp_base_opt(:,3)'==2);
    posError_T_minor = find(align_vec_mat(i,:)==0 & tmp_base_opt(:,4)'==2);
    
    % the Sanger sequence errors
    
    % first line: amplified sequence at position
    writeErrors(1,[posError_A_major,posError_C_major,posError_G_major,posError_T_major]) = [s1(posError_A_major),s1(posError_C_major),s1(posError_G_major),s1(posError_T_major)]; % amplified
    writeErrors(1,[posError_A_minor,posError_C_minor,posError_G_minor,posError_T_minor]) = [s1(posError_A_minor),s1(posError_C_minor),s1(posError_G_minor),s1(posError_T_minor)]; % amplified
                                                                                                                                  
    % second line - Sanger major allele
    writeErrors(2,[posError_A_major,posError_C_major,posError_G_major,posError_T_major]) =  [s2(posError_A_major),s2(posError_C_major),s2(posError_G_major),s2(posError_T_major)];% major allele at errors

    % third line - Sanger minor allele
    writeErrors(3,posError_A_minor) = 'A';
    writeErrors(3,posError_C_minor) = 'C';
    writeErrors(3,posError_G_minor) = 'G';
    writeErrors(3,posError_T_minor) = 'T';
    
    
    s1_to_match.writeErrors{i} = writeErrors;
    
    
    
    
    
    % align sequences of 750 and 450
    for j=1:length(res{i})
      
      % 750
      if ~isempty(res{i}{j}.s750.Header)
        
        % 750
        % collect information in gloabal coordinates of the full
        len_curr_full = length(res{i}{j}.Full.Sequence);
        coor = char('P'*ones(3,len_curr_full));
        coor(1,res{i}{j}.s750.Start(1):res{i}{j}.s750.Stop(1)) = res{i}{j}.s750.Sequence;
        start_pos = res{i}{j}.Full.Start(2)+AlignStruct.Start(1)-1+(AmplifiedSeq_start-1);
        pos_s1_s2 = start_pos:start_pos+curr_lengthValidSanger;
        coor(2,pos_s1_s2) = s1;
        % part of Sanger to be compared to amplified
        coor(3,pos_s1_s2) = s2;
        
        tmpResMajor = 10*ones(1,len_curr_full);
        tmpResMajor(pos_s1_s2) = align_mat(i,:);
        
        tmpResMinor = 10*ones(1,len_curr_full);
        tmpResMinor(pos_s1_s2) = align_vec_mat(i,:);
        
        % length of overlap between Sanger and 750 
        % number of errors in that region - with the major and the minor
        overlapIndices = find(coor(1,:)~='P' & coor(3,:)~='P');
        overlapData{i}{j}.s750.length = length(overlapIndices);
        overlapData{i}{j}.s750.numErrorsMajor = length(find(tmpResMajor(overlapIndices)==0)); 
        overlapData{i}{j}.s750.numErrorsMinor = length(find(tmpResMinor(overlapIndices)<1)); 
      
      
      end % 750
      
      % 450
      if ~isempty(res{i}{j}.s450.Header)
        
        % 450
        % collect information in gloabal coordinates of the full
        len_curr_full = length(res{i}{j}.Full.Sequence);
        coor = char('P'*ones(3,len_curr_full));
        coor(1,res{i}{j}.s450.Start(1):res{i}{j}.s450.Stop(1)) = res{i}{j}.s450.Sequence;
        start_pos = res{i}{j}.Full.Start(2)+AlignStruct.Start(1)-1+(AmplifiedSeq_start-1);
        pos_s1_s2 = start_pos:start_pos+curr_lengthValidSanger;
        coor(2,pos_s1_s2) = s1;
        % part of Sanger to be compared to amplified
        coor(3,pos_s1_s2) = s2;
        
        tmpResMajor = 10*ones(1,len_curr_full);
        tmpResMajor(pos_s1_s2) = align_mat(i,:);
        
        tmpResMinor = 10*ones(1,len_curr_full);
        tmpResMinor(pos_s1_s2) = align_vec_mat(i,:);
        
        % length of overlap between Sanger and 450 
        % number of errors in that region - with the major and the minor
        overlapIndices = find(coor(1,:)~='P' & coor(3,:)~='P');
        overlapData{i}{j}.s450.length = length(overlapIndices);
        overlapData{i}{j}.s450.numErrorsMajor = length(find(tmpResMajor(overlapIndices)==0)); 
        overlapData{i}{j}.s450.numErrorsMinor = length(find(tmpResMinor(overlapIndices)<1)); 
      
      end
      
      
    end
    
    
  end
else % Flag==2
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %keyboard
  for i=1:length(res)
    %keyboard
        
    % amplified
    s1 = seqrcomplement(res{i}{1}.amp.Sequence(max(1,length(res{i}{1}.amp.Sequence)-len_seqSanger):end-(AmplifiedSeq_start-1)));
    s2 = seqSanger(SangerSeq_start:len_seqSanger); % Sanger major allele 
    
    AlignStruct = localalign(s1, s2,'Alphabet','NT');
    AlignStruct.Readme = 'the first is amplified';
    
    
    s1 = s1(AlignStruct.Start(1):(AlignStruct.Start(1)+curr_lengthValidSanger)); % s1 is the amplified to later compared to
    s2 = s2(AlignStruct.Start(2):(AlignStruct.Start(2)+curr_lengthValidSanger));
    
    %keyboard
    tmp_base_opt = base_opt(SangerSeq_start:len_seqSanger,:);
    tmp_base_opt = tmp_base_opt(AlignStruct.Start(2):(AlignStruct.Start(2)+curr_lengthValidSanger),:);
    [dist(i),align_vec_mat(i,:)] = seqDistanceWithAmbiguity(s1,tmp_base_opt);

    %keyboard
    
    % align_vec_mat is 1 when sequence is similar to Sanger major allele or Sanger minor allele
    align_mat(i,:) = (s1==s2); % matrix of 1 when sequence matches Sanger major allele and 0 otherwise
    align_vec_mat(i,align_vec_mat(i,:)>align_mat(i,:)) = 2; % change minor to major
    % matched amplified sequence to later compare
    s1_to_match.seq{i} = s1;
    s1_to_match.SangerOffset(i) = SangerSeq_start+AlignStruct.Start(2)-1;
    
    writeErrors =  char('P'*ones(3,length(s1)));
    
    posError_A_major = find(align_vec_mat(i,:)==0 & tmp_base_opt(:,1)'==1);
    posError_C_major = find(align_vec_mat(i,:)==0 & tmp_base_opt(:,2)'==1);
    posError_G_major = find(align_vec_mat(i,:)==0 & tmp_base_opt(:,3)'==1);
    posError_T_major = find(align_vec_mat(i,:)==0 & tmp_base_opt(:,4)'==1);
    
    posError_A_minor = find(align_vec_mat(i,:)==0 & tmp_base_opt(:,1)'==2);
    posError_C_minor = find(align_vec_mat(i,:)==0 & tmp_base_opt(:,2)'==2);
    posError_G_minor = find(align_vec_mat(i,:)==0 & tmp_base_opt(:,3)'==2);
    posError_T_minor = find(align_vec_mat(i,:)==0 & tmp_base_opt(:,4)'==2);
    
    
    
    
    % the Sanger sequence errors
    
    % first line: amplified sequence at position
    writeErrors(1,[posError_A_major,posError_C_major,posError_G_major,posError_T_major]) = [s1(posError_A_major),s1(posError_C_major),s1(posError_G_major),s1(posError_T_major)]; % amplified
    writeErrors(1,[posError_A_minor,posError_C_minor,posError_G_minor,posError_T_minor]) = [s1(posError_A_minor),s1(posError_C_minor),s1(posError_G_minor),s1(posError_T_minor)]; % amplified
                                                                                                                                  
    % second line - Sanger major allele
    writeErrors(2,[posError_A_major,posError_C_major,posError_G_major,posError_T_major]) =  [s2(posError_A_major),s2(posError_C_major),s2(posError_G_major),s2(posError_T_major)];% major allele at errors

    % third line - Sanger minor allele
    writeErrors(3,posError_A_minor) = 'A';
    writeErrors(3,posError_C_minor) = 'C';
    writeErrors(3,posError_G_minor) = 'G';
    writeErrors(3,posError_T_minor) = 'T';
    
    
    s1_to_match.writeErrors{i} = writeErrors;
    
   
    
    % align sequences of 750 and 450
    for j=1:length(res{i})
      
      % 750
      if ~isempty(res{i}{j}.s750.Header)
        
        % collect information in gloabal coordinates of the full
        len_curr_full = length(res{i}{j}.Full.Sequence);
        coor_revcomp = char('P'*ones(3,len_curr_full));
        % 750
        coor_revcomp(1,len_curr_full-res{i}{j}.s750.Stop(1):len_curr_full-res{i}{j}.s750.Start(1)) = seqrcomplement(res{i}{j}.s750.Sequence);
        
        % path of amplified to compare to Sanger
        start_pos = (len_curr_full-res{i}{j}.Full.Stop(2))+AlignStruct.Start(1)-1+(AmplifiedSeq_start-1);
        pos_s1_s2 = [start_pos:start_pos+curr_lengthValidSanger];
        coor_revcomp(2,pos_s1_s2) = s1;
        % part of Sanger to be compared to amplified
        coor_revcomp(3,pos_s1_s2) = s2;
        
        
        tmpResMajor = 10*ones(1,len_curr_full);
        tmpResMajor(pos_s1_s2) = align_mat(i,:);
        
        tmpResMinor = 10*ones(1,len_curr_full);
        tmpResMinor(pos_s1_s2) = align_vec_mat(i,:);
        
        % length of overlap between Sanger and 750 
        % number of errors in that region - with the major and the minor
        overlapIndices = find(coor_revcomp(1,:)~='P' & coor_revcomp(3,:)~='P');
        overlapData{i}{j}.s750.length = length(overlapIndices);
        overlapData{i}{j}.s750.numErrorsMajor = length(find(tmpResMajor(overlapIndices)==0)); 
        overlapData{i}{j}.s750.numErrorsMinor = length(find(tmpResMinor(overlapIndices)<1)); 
      end % 750
      
      % 450
      if ~isempty(res{i}{j}.s450.Header)
        % collect information in gloabal coordinates of the full
        len_curr_full = length(res{i}{j}.Full.Sequence);
        coor_revcomp = char('P'*ones(3,len_curr_full));
        % 450
        coor_revcomp(1,len_curr_full-res{i}{j}.s450.Stop(1):len_curr_full-res{i}{j}.s450.Start(1)) = seqrcomplement(res{i}{j}.s450.Sequence);
        
        % path of amplified to compare to Sanger
        start_pos = (len_curr_full-res{i}{j}.Full.Stop(2))+AlignStruct.Start(1)-1+(AmplifiedSeq_start-1);
        pos_s1_s2 = [start_pos:start_pos+curr_lengthValidSanger];
        coor_revcomp(2,pos_s1_s2) = s1;
        % part of Sanger to be compared to amplified
        coor_revcomp(3,pos_s1_s2) = s2;
        
        
        tmpResMajor = 10*ones(1,len_curr_full);
        tmpResMajor(pos_s1_s2) = align_mat(i,:);
        
        tmpResMinor = 10*ones(1,len_curr_full);
        tmpResMinor(pos_s1_s2) = align_vec_mat(i,:);
        
        % length of overlap between Sanger and 450 
        % number of errors in that region - with the major and the minor
        overlapIndices = find(coor_revcomp(1,:)~='P' & coor_revcomp(3,:)~='P');
        overlapData{i}{j}.s450.length = length(overlapIndices);
        overlapData{i}{j}.s450.numErrorsMajor = length(find(tmpResMajor(overlapIndices)==0)); 
        overlapData{i}{j}.s450.numErrorsMinor = length(find(tmpResMinor(overlapIndices)<1)); 
      end
      
      
    end
    
    
  end
  
  
  
end
%disp('is start_s1=1 ok?');pause

% extract

% align_mat=align_vec_mat=0 mean equal to major
% align_vec_mat=1 & align_mat=0 - equal to minor but not to major in the sequence
% align_mat=1 & align_vec_mat==0 can not be

if ~isempty(find(align_mat==1 & align_vec_mat==0))
  disp('problem');keyboard
end

% find amplified that match the sanger aligned regions
distMatOfComparedSeqs = squareform(seqpdist(s1_to_match.seq,'Alphabet','NT','Method','p-distance','PairwiseAlignment',false));
[xx,yy] = find(triu(distMatOfComparedSeqs,1)==0);
del = find(xx>=yy);
xx(del) = [];
yy(del) = [];

classSeqs = unionfind([xx,yy],size(distMatOfComparedSeqs,1));

[vals inds num_dups] = get_duplicates(classSeqs);
lines = zeros(length(inds),1);
%keyboard
SangerOffset = cell(length(inds),1);
writeErrors = cell(length(inds),1);
for i=1:length(inds)
  lines(i) = inds{i}(1);
  SangerOffset{i} = s1_to_match.SangerOffset(inds{i});
  writeErrors{i} = s1_to_match.writeErrors(inds{i});
end

dist = dist(lines);
align_vec_mat = align_vec_mat(lines,:);
align_mat = align_mat(lines,:);

%keyboard

[~,xi] = sort(dist);

% use one row from each class
align_vec_mat = align_vec_mat(xi,:); % sort the rows
align_mat = align_mat(xi,:);
inds = inds(xi);
SangerOffset = SangerOffset(xi);
writeErrors = writeErrors(xi);

%keyboard

% color code: 
% black - (0) - sequence different from major allele and of minor allele (if exists)
% gray - (1) - sequence same as Sanger major allele
% white - (2) - sequence same as minor allele








