%split 750 450  - what is the order we make? open names

function results_chroma3_in(basicDirName,fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,tmp_ind_W_454,tmp_ind_W_Solexa,sameNumbers,paperPosition,thresh,fs,th_sec2first)
%keyboar
lenName = 50;

if ~isempty(findstr(fileName,'F.scf'))
  F_flag = 1;
else
  F_flag = 2; % reverse
end
%keyboard

% find Sanger sequences

[seqSanger,base_opt] = find_possible_seq_from_chromatogram_noam3([basicDirName,fileName],th_sec2first);

% SeqS = seqSanoger;SB = base_opt;save ~/tmp SeqS SB
%keyboard
% collect data and errors
[align_vec_mat,align_mat,overlapData,res,inds,SangerOffset,writeErrors] = compare_sangerseq_to_ampseq2(res,lengthValidSanger,seqSanger,F_flag,base_opt,SangerSeq_start,AmplifiedSeq_start);


%keyboard
collect_number750 = [];
collect_number450 = [];
for i=1:length(inds)
  %if i==4,keyboard,end
  ampData = [];
  lab750 = [];
  lab450 = [];
  for j=1:length(inds{i})  
    
    %if res{inds{i}(j)}{1}.Full.Numbers==220048,keyboard,end 
    
    % amplified name and number of errors - with minor
    ampData{j} = [];
    for k=1:length(res{inds{i}(j)})
      for l=1:length(res{inds{i}(j)}{k}.Full.Numbers)
        ampData{j} = [ampData{j},' {\color{magenta}',num2str(res{inds{i}(j)}{k}.Full.Numbers(l)),'} '];
      end
      
    end
    ampData{j} = [ampData{j}, ' {\color{magenta} [',num2str(size(align_mat,2)),' ',num2str([length(find(align_mat(i,:)==0)) length(find(align_vec_mat(i,:)==0))]),']}'];
    
    ampData{j} = delSpaces(ampData{j});
    
    
    % 750
    clear data750 number750
    for k=1:length(res{inds{i}(j)})
      if ~isempty(res{inds{i}(j)}{k}.s750.Sequence)
        if ~isempty(intersect(res{inds{i}(j)}{k}.s750.Numbers,tmp_ind_W_Solexa))
%keyboard
          number750{k} = res{inds{i}(j)}{k}.s750.Numbers;
          data750{k} = [ overlapData{inds{i}(j)}{k}.s750.length overlapData{inds{i}(j)}{k}.s750.numErrorsMajor overlapData{inds{i}(j)}{k}.s750.numErrorsMinor]; 
        else
          number750{k} = [];
          data750{k} = [];
        end
      else
        number750{k} = [];
        data750{k} = [];
      end
    end
    
    taken750 = [];
    tmp_lab750 = [];
    for k=1:length(res{inds{i}(j)})
      if ~isempty(number750{k}) && isempty(intersect(taken750,number750{k}))
        for l=1:length(number750{k})
          tmp_lab750 = [tmp_lab750,';{\color{blue} {\bf ',num2str(number750{k}(l)'),'}}'];
        end
        
        tmp_lab750 = [tmp_lab750,';{\color{blue} {\bf [',num2str(data750{k}),']}}'];

        taken750 = [taken750,number750{k}'];
      end
    end
    
    lab750{j} = delSpaces(tmp_lab750);
    collect_number750 = [collect_number750,taken750];
    
    % 450
    clear data450 number450
    for k=1:length(res{inds{i}(j)})
      if ~isempty(res{inds{i}(j)}{k}.s450.Sequence)
        if ~isempty(intersect(res{inds{i}(j)}{k}.s450.Numbers,tmp_ind_W_454))
%keyboard
          number450{k} = res{inds{i}(j)}{k}.s450.Numbers;
          data450{k} = [ overlapData{inds{i}(j)}{k}.s450.length overlapData{inds{i}(j)}{k}.s450.numErrorsMajor overlapData{inds{i}(j)}{k}.s450.numErrorsMinor]; 
        else
          number450{k} = [];
          data450{k} = [];
        end
      else
        number450{k} = [];
        data450{k} = [];
      end
    end
    
    taken450 = [];
    tmp_lab450 = [];
    for k=1:length(res{inds{i}(j)})
      if ~isempty(number450{k}) && isempty(intersect(taken450,number450{k}))
        for l=1:length(number450{k})
          tmp_lab450 = [tmp_lab450,';{\color{red} {\bf ',num2str(number450{k}(l)'),'}}'];
        end
        
        tmp_lab450 = [tmp_lab450,';{\color{blue} {\bf [',num2str(data450{k}),']}}'];

        taken450 = [taken450,number450{k}'];
      end
    end
    
    lab450{j} = delSpaces(tmp_lab450);
    collect_number450 = [collect_number450,taken450];
    
    
    
    
    
  end
    
  
 
   
  
  lab{i} = [];
  for j=1:length(inds{i})  
    lab{i} = [lab{i},';',ampData{j}];
    lab{i} = [lab{i},';',lab750{j},';',lab450{j}];
  end
  
  
  %keyboard
  % write the text for the errors
  for j=1:length(inds{i})
    findErrors = find(writeErrors{i}{j}(1,:)~='P');
    textErrors{i}.position = findErrors;
    
    for k=1:length(findErrors)
      currPos = findErrors(k)+SangerOffset{i}(j)-1;
      if currPos>=100
        positionInText = num2str(currPos)';
      elseif currPos>=10 & currPos<100
        positionInText = [0;num2str(currPos)'];
      else % <10
        positionInText = [0;0;num2str(currPos)'];
      end
      textErrors{i}.valAtError(:,k) = [writeErrors{i}{j}(:,findErrors(k));positionInText];
    end
  end
  %keyboard
  
end


%keyboard
for i=1:length(lab)
  if length(lab{i})>thresh
    tmpLab = lab{i};
    newLab = [];
    while ~isempty(tmpLab)
      a = findstr(tmpLab,';{\color');
      b = find(a>thresh);
      if ~isempty(b)
        newLab = [newLab,tmpLab(1:a(b(1))-1),' \newline '];
        tmpLab(1:a(b(1))-1) = [];
      else
        newLab = [newLab,tmpLab];
        tmpLab = [];
      end
      
    end

    lab{i} = newLab;
  end
end

%keyboard
disp('compare colors to former one!!!')
%keyboard
% find those that were not amplified and appear in our list
collect_number750 = unique(collect_number750);
not_amplified_Solexa = setdiff(tmp_ind_W_Solexa,collect_number750);
not_amplified_454 = setdiff(tmp_ind_W_454,collect_number450);


%%%%%%%%
% plot



% add those that were not amplified


align_vec_mat = [align_vec_mat;2*ones(length(not_amplified_Solexa),size(align_vec_mat,2))];
add_yl_Solexa = cell(1,length(not_amplified_Solexa));
for i=1:length(not_amplified_Solexa)
  add_yl_Solexa{i} = ['{\color{blue} ',num2str(not_amplified_Solexa(i)),' NA}'];
end

align_vec_mat = [align_vec_mat;2*ones(length(not_amplified_454),size(align_vec_mat,2))];
add_yl_454 = cell(1,length(not_amplified_454));
for i=1:length(not_amplified_454)
  add_yl_454{i} = ['{\color{red} ',num2str(not_amplified_454(i)),' NA}'];
end

xl = cell(1,length(pos_xl))
for i=1:length(xl)
  xl{i} = num2str(pos_xl(i));
end
%keyboard

clf
imagesc(align_vec_mat);
for i=1:size(align_vec_mat,1)
  line([1 size(align_vec_mat,2)],[i+0.5,i+0.5]','color','k');        
end

for i=1:length(textErrors)
  for j=1:length(textErrors{i}.position)
    text(textErrors{i}.position(j),i,textErrors{i}.valAtError(:,j),'color','r','fontsize',8);
  end
end


colormap('gray')    

if ~exist('fs')
  fs = 10;
end

toPrint = [lab,add_yl_454,add_yl_Solexa];
format_ticks(gca,xl,toPrint,[pos_xl],[1:size(align_vec_mat,1)],[],[],[],'fontweight','bold','fontsize',fs)
name = fileName;
name = strrep(name,'_','\_')
title(name,'fontweight','bold','fontsize',fs+10)
%keyboard

set(gcf,'paperorientation','landscape','paperposition',[-1    0    12    8]); %[from left from_bottom XX height ]
set(gca,'fontweight','bold','fontsize',fs,'position',paperPosition)

%keyboard