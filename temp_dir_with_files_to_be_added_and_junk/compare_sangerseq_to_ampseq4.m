function [align_mat,align_vec_mat,overlapData,saveErrors,SangerOffset,reordered_cData,errorsMinorOverSequence]=compare_sangerseq_to_ampseq4(cData,lengthValidSanger,seqSanger,F_flag,base_opt,SangerSeq_start,AmplifiedSeq_start)

len_seqSanger = length(seqSanger);

curr_lengthValidSanger = min([lengthValidSanger,len_seqSanger-AmplifiedSeq_start-2]);

% n=400;
align_mat=zeros(length(cData),curr_lengthValidSanger+1);
align_vec_mat=zeros(length(cData),curr_lengthValidSanger+1);
dist=zeros(length(cData),1);

overlapData = struct;
overlapData.amp = zeros(length(cData.amp),3);
overlapData.s750 = overlapData.amp;
overlapData.s450 = overlapData.amp;

SangerOffset = zeros(length(cData.amp),1);

%keyboard
if F_flag==1
  for i=1:length(cData.amp)
   
    %if cData.Numbers(i)==579340
    %  keyboard
    %end
    
    clear tmp_align_vec_mat tmp_align_mat tmp_dist numMinorErrors store_base_opt tmp_AlignStruct tmp_s1 tmp_s2 store_s1 store_s2
    for j=1:length(cData.amp{i}.Sequence)
       
      % amplified
       tmp_s1 =  cData.amp{i}.Sequence{j}(AmplifiedSeq_start:len_seqSanger);
       tmp_s2 = seqSanger(SangerSeq_start:len_seqSanger); % Sanger major allele 
       
       tmp_AlignStruct{j} = localalign(tmp_s1, tmp_s2,'Alphabet','NT');
       tmp_AlignStruct{j}.Readme = 'the first is amplified';
       
       tmp_s1 = tmp_s1(tmp_AlignStruct{j}.Start(1):(tmp_AlignStruct{j}.Start(1)+curr_lengthValidSanger)); % s1 is the amplified to later compared to
       tmp_s2 = tmp_s2(tmp_AlignStruct{j}.Start(2):(tmp_AlignStruct{j}.Start(2)+curr_lengthValidSanger));
       
       store_base_opt{j} = base_opt(SangerSeq_start:len_seqSanger,:);
       store_base_opt{j} = store_base_opt{j}(tmp_AlignStruct{j}.Start(2):(tmp_AlignStruct{j}.Start(2)+curr_lengthValidSanger),:);
       
       store_s1{j} = tmp_s1;
       store_s2{j} = tmp_s2;
       
       
       [tmp_dist(i),tmp_align_vec_mat(j,:)] = seqDistanceWithAmbiguity(tmp_s1,store_base_opt{j});
       % align_vec_mat is 1 when sequence is similar to Sanger major allele or Sanger minor allele
       tmp_align_mat(j,:) = (tmp_s1==tmp_s2);% matrix of 1 when sequence matches Sanger major allele and 0 otherwise
       
       tmp_align_vec_mat(j,tmp_align_vec_mat(j,:)>tmp_align_mat(j,:)) = 2;% change minor to major
       
       numMinorErrors(j) = length(find(tmp_align_vec_mat(j,:)==0));
       
    end
    
    % the the amplified which best matches
    if i==1
     disp('taking the amplified which best matches')
    end
    
    [junk,indMinMinorError] = min(numMinorErrors);
    align_mat(i,:) = tmp_align_mat(indMinMinorError,:);
    align_vec_mat(i,:) = tmp_align_vec_mat(indMinMinorError,:);
    AlignStruct = tmp_AlignStruct{indMinMinorError};
    s1 = store_s1{indMinMinorError};
    s2 = store_s2{indMinMinorError};
    tmp_base_opt = store_base_opt{indMinMinorError};
    
    
    
    minorAlleleErrorPosition = find(align_vec_mat(i,:)==0);
    
    overlapData.amp(i,:) = [size(align_vec_mat,2) length(find(align_mat(i,:)==0)) length(minorAlleleErrorPosition)];
    
    %%%%%%%%%%%%%%%
    % Errors
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
    
    saveErrors{i} = writeErrors; 
    
    % Sanger offset
    SangerOffset(i) = SangerSeq_start+AlignStruct.Start(2)-1;
    
    %%%%%%%%%%%%% save data of amplified sequence - state where the minor allele was wrong
    len_curr_full = length(cData.Full{i}.Sequence);
    coor = char('P'*ones(2,len_curr_full));
    coor(1,:) = cData.Full{i}.Sequence;
    start_pos = cData.amp{i}.Start(indMinMinorError,2)+AlignStruct.Start(1)-1+(AmplifiedSeq_start-1);
    pos_s1_s2 = start_pos:start_pos+curr_lengthValidSanger;
    % part of amplified to compare - enter the minor allele
    currAllele = s1;
    currAllele(minorAlleleErrorPosition) = 'N';
    coor(2,pos_s1_s2) = currAllele;
    errorsMinorOverSequence{i} = coor;clear coor
    % end save data of amplified sequence
    if i==1
     disp('print N where minor allele is incorrect')
    end
    
    % align sequences of 750 and 450
    if ~isempty(cData.s750{i}.Header)
      
      % 750
      % collect information in global coordinates of the full
      len_curr_full = length(cData.Full{i}.Sequence);
      coor = char('P'*ones(3,len_curr_full));
      coor(1,cData.s750{i}.Start(2):cData.s750{i}.Stop(2)) = cData.s750{i}.Sequence;
      start_pos = cData.amp{i}.Start(indMinMinorError,2)+AlignStruct.Start(1)-1+(AmplifiedSeq_start-1);
      pos_s1_s2 = start_pos:start_pos+curr_lengthValidSanger;
      coor(2,pos_s1_s2) = s1;
      % part of Sanger to be compared to amplified
      coor(3,pos_s1_s2) = s2;
      
      tmpcDataMajor = 10*ones(1,len_curr_full);
      tmpcDataMajor(pos_s1_s2) = align_mat(i,:);
      
      tmpcDataMinor = 10*ones(1,len_curr_full);
      tmpcDataMinor(pos_s1_s2) = align_vec_mat(i,:);
      
      % length of overlap between Sanger and 750 
      % number of errors in that region - with the major and the minor
      overlapIndices = find(coor(1,:)~='P' & coor(3,:)~='P');
      overlapData.s750(i,:) = [length(overlapIndices) length(find(tmpcDataMajor(overlapIndices)==0)) length(find(tmpcDataMinor(overlapIndices)<1))];
    else
      overlapData.s750(i,:) = [-1 -1 -1];
    end % 750
    
    
    
    if ~isempty(cData.s450{i}.Header)
      
      % 450
      % collect information in gloabal coordinates of the full
      len_curr_full = length(cData.Full{i}.Sequence);
      coor = char('P'*ones(3,len_curr_full));
      coor(1,cData.s450{i}.Start(2):cData.s450{i}.Stop(2)) = cData.s450{i}.Sequence;
      start_pos = cData.amp{i}.Start(indMinMinorError,2)+AlignStruct.Start(1)-1+(AmplifiedSeq_start-1);
      pos_s1_s2 = start_pos:start_pos+curr_lengthValidSanger;
      coor(2,pos_s1_s2) = s1;
      % part of Sanger to be compared to amplified
      coor(3,pos_s1_s2) = s2;
      
      tmpcDataMajor = 10*ones(1,len_curr_full);
      tmpcDataMajor(pos_s1_s2) = align_mat(i,:);
      
      tmpcDataMinor = 10*ones(1,len_curr_full);
      tmpcDataMinor(pos_s1_s2) = align_vec_mat(i,:);
      
      % length of overlap between Sanger and 450 
      % number of errors in that region - with the major and the minor
      overlapIndices = find(coor(1,:)~='P' & coor(3,:)~='P');
      overlapData.s450(i,:) = [length(overlapIndices) length(find(tmpcDataMajor(overlapIndices)==0)) length(find(tmpcDataMinor(overlapIndices)<1))];
    else
      overlapData.s450(i,:) = [-1 -1 -1];
    end % 450
    
    
    
  end

else % Flag==2
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %keyboard
  
  for i=1:length(cData.amp)
     
    
    
    clear tmp_align_vec_mat tmp_align_mat tmp_dist numMinorErrors store_base_opt tmp_AlignStruct tmp_s1 tmp_s2 store_s1 store_s2
    for j=1:length(cData.amp{i}.Sequence)
      
      % amplified
       tmp_s1 = seqrcomplement(cData.amp{i}.Sequence{j}(max(1,length(cData.amp{i}.Sequence{j})-len_seqSanger):end-(AmplifiedSeq_start-1)));
       tmp_s2 = seqSanger(SangerSeq_start:len_seqSanger); % Sanger major allele 
       
       tmp_AlignStruct{j} = localalign(tmp_s1, tmp_s2,'Alphabet','NT');
       tmp_AlignStruct{j}.Readme = 'the first is amplified';
       
       tmp_s1 = tmp_s1(tmp_AlignStruct{j}.Start(1):(tmp_AlignStruct{j}.Start(1)+curr_lengthValidSanger)); % s1 is the amplified to later compared to
       tmp_s2 = tmp_s2(tmp_AlignStruct{j}.Start(2):(tmp_AlignStruct{j}.Start(2)+curr_lengthValidSanger));
       
       store_base_opt{j} = base_opt(SangerSeq_start:len_seqSanger,:);
       store_base_opt{j} = store_base_opt{j}(tmp_AlignStruct{j}.Start(2):(tmp_AlignStruct{j}.Start(2)+curr_lengthValidSanger),:);
       
       store_s1{j} = tmp_s1;
       store_s2{j} = tmp_s2;
       
       [tmp_dist(i),tmp_align_vec_mat(j,:)] = seqDistanceWithAmbiguity(tmp_s1,store_base_opt{j});
       % align_vec_mat is 1 when sequence is similar to Sanger major allele or Sanger minor allele
       tmp_align_mat(j,:) = (tmp_s1==tmp_s2);% matrix of 1 when sequence matches Sanger major allele and 0 otherwise
       tmp_align_vec_mat(j,tmp_align_vec_mat(j,:)>tmp_align_mat(j,:)) = 2;% change minor to major       
       numMinorErrors(j) = length(find(tmp_align_vec_mat(j,:)==0));
    end
    
    % the the amplified which best matches
    if i==1
     disp('taking the amplified which best matches')
    end
    
    [junk,indMinMinorError] = min(numMinorErrors);
    align_mat(i,:) = tmp_align_mat(indMinMinorError,:);
    align_vec_mat(i,:) = tmp_align_vec_mat(indMinMinorError,:);
    AlignStruct = tmp_AlignStruct{indMinMinorError};
    s1 = store_s1{indMinMinorError};
    s2 = store_s2{indMinMinorError};
    tmp_base_opt = store_base_opt{indMinMinorError};
    
    minorAlleleErrorPosition = find(align_vec_mat(i,:)==0);
    
    overlapData.amp(i,:) = [size(align_vec_mat,2) length(find(align_mat(i,:)==0)) length(minorAlleleErrorPosition)];
    
    
    % Errors
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
    
    saveErrors{i} = writeErrors; 
        
    % Sanger offset
    SangerOffset(i) = SangerSeq_start+AlignStruct.Start(2)-1;
    %keyboard
    
    
    %%%%%%%%%%%%% save data of amplified sequence - state where the minor allele was wrong
    len_curr_full = length(cData.Full{i}.Sequence);
    coor_revcomp = char('P'*ones(2,len_curr_full));
    % 750
    coor_revcomp(1,:) = seqrcomplement(cData.Full{i}.Sequence);
    % path of amplified to compare to Sanger
    start_pos = (len_curr_full-cData.amp{i}.Stop(indMinMinorError,2))+AlignStruct.Start(1)+(AmplifiedSeq_start-1);
    pos_s1_s2 = [start_pos:start_pos+curr_lengthValidSanger];

    % part of amplified to compare - enter the minor allele
    currAllele = s1;
    currAllele(minorAlleleErrorPosition) = 'N';
    coor_revcomp(2,pos_s1_s2) = currAllele;
    if i==1
     disp('print N where minor allele is incorrect. Assumes no non-ACGT')
    end
    
    %keyboard
    clear coor
    coor(1,:) = seqrcomplement(coor_revcomp(1,:));
    coor(1,:) = strrep(coor(1,:),'*','P');
    coor(2,:) = seqrcomplement(coor_revcomp(2,:));
    coor(2,:) = strrep(coor(2,:),'*','P');
    errorsMinorOverSequence{i} = coor;clear coor

    % end save data of amplified sequence
    % align sequences of 750 and 450
    if ~isempty(cData.s750{i}.Header)
      % collect information in gloabal coordinates of the full
       len_curr_full = length(cData.Full{i}.Sequence);
       coor_revcomp = char('P'*ones(3,len_curr_full));
       % 750
       coor_revcomp(1,len_curr_full-cData.s750{i}.Stop(2):len_curr_full-cData.s750{i}.Start(2)) = seqrcomplement(cData.s750{i}.Sequence);
       
       % path of amplified to compare to Sanger
       start_pos = (len_curr_full-cData.amp{i}.Stop(indMinMinorError,2))+AlignStruct.Start(1)-1+(AmplifiedSeq_start-1);
       pos_s1_s2 = [start_pos:start_pos+curr_lengthValidSanger];
       coor_revcomp(2,pos_s1_s2) = s1;
       % part of Sanger to be compared to amplified
       coor_revcomp(3,pos_s1_s2) = s2;
       
       tmpcDataMajor = 10*ones(1,len_curr_full);
       tmpcDataMajor(pos_s1_s2) = align_mat(i,:);
       
       tmpcDataMinor = 10*ones(1,len_curr_full);
       tmpcDataMinor(pos_s1_s2) = align_vec_mat(i,:);
       
       % length of overlap between Sanger and 750 
       % number of errors in that region - with the major and the minor
       overlapIndices = find(coor_revcomp(1,:)~='P' & coor_revcomp(3,:)~='P');
       overlapData.s750(i,:) = [length(overlapIndices) length(find(tmpcDataMajor(overlapIndices)==0)) length(find(tmpcDataMinor(overlapIndices)<1))];
    else
      overlapData.s750(i,:) = [-1 -1 -1];
    end % 750
    
    
    
   % 450
   if ~isempty(cData.s450{i}.Header) 
     % collect information in gloabal coordinates of the full
       len_curr_full = length(cData.Full{i}.Sequence);
       coor_revcomp = char('P'*ones(3,len_curr_full));
       % 450
       coor_revcomp(1,len_curr_full-cData.s450{i}.Stop(2):len_curr_full-cData.s450{i}.Start(2)) = seqrcomplement(cData.s450{i}.Sequence);
       
       % path of amplified to compare to Sanger
       start_pos = (len_curr_full-cData.amp{i}.Stop(indMinMinorError,2))+AlignStruct.Start(1)-1+(AmplifiedSeq_start-1);
       pos_s1_s2 = [start_pos:start_pos+curr_lengthValidSanger];
       coor_revcomp(2,pos_s1_s2) = s1;
       % part of Sanger to be compared to amplified
       coor_revcomp(3,pos_s1_s2) = s2;
       
       tmpcDataMajor = 10*ones(1,len_curr_full);
       tmpcDataMajor(pos_s1_s2) = align_mat(i,:);
       
       tmpcDataMinor = 10*ones(1,len_curr_full);
       tmpcDataMinor(pos_s1_s2) = align_vec_mat(i,:);
       
       % length of overlap between Sanger and 450 
       % number of errors in that region - with the major and the minor
       overlapIndices = find(coor_revcomp(1,:)~='P' & coor_revcomp(3,:)~='P');
       overlapData.s450(i,:) = [length(overlapIndices) length(find(tmpcDataMajor(overlapIndices)==0)) length(find(tmpcDataMinor(overlapIndices)<1))];
   else
      overlapData.s450(i,:) = [-1 -1 -1];
   end
  end
  
end

% align_mat=align_vec_mat=0 mean equal to major
% align_vec_mat=1 & align_mat=0 - equal to minor but not to major in the sequence
% align_mat=1 & align_vec_mat==0 can not be

if ~isempty(find(align_mat==1 & align_vec_mat==0))
  disp('problem');keyboard
end

[junk,indSort] = sortrows(overlapData.amp,[3,2]);

overlapData.amp = overlapData.amp(indSort,:);
overlapData.s750 = overlapData.s750(indSort,:);
overlapData.s450 = overlapData.s450(indSort,:);
%keyboard


align_mat = align_mat(indSort,:);
align_vec_mat = align_vec_mat(indSort,:);

saveErrors = saveErrors(indSort);
SangerOffset = SangerOffset(indSort);

errorsMinorOverSequence = errorsMinorOverSequence(indSort);

reordered_cData = struct;
reordered_cData.Numbers = cData.Numbers(indSort);
reordered_cData.Full = cData.Full(indSort);
reordered_cData.amp = cData.amp(indSort);
reordered_cData.s750 = cData.s750(indSort);
reordered_cData.s450 = cData.s450(indSort);


%keyboard

% color code: 
% black - (0) - sequence different from major allele and of minor allele (if exists)
% gray - (1) - sequence same as Sanger major allele
% white - (2) - sequence same as minor allele








