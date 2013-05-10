%split 750 450  - what is the order we make? open names

function resData=results_chroma2_in(basicDirName,fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,tmp_ind_W_454,tmp_ind_W_Solexa,sameNumbers,paperPosition,thresh,fs)
%keyboar
lenName = 50;

if ~isempty(findstr(fileName,'F.scf'))
  F_flag = 1;
else
  F_flag = 2; % reverse
end
%keyboard

% find Sanger sequences
[seqSanger,base_opt] = find_possible_seq_from_chromatogram_noam([basicDirName,fileName]);

% SeqS = seqSanoger;SB = base_opt;save ~/tmp SeqS SB
%keyboard
% collect data and errors
[align_vec_mat,align_mat,overlapData,res,inds,SangerOffset,writeErrors] = compare_sangerseq_to_ampseq2(res,lengthValidSanger,seqSanger,F_flag,base_opt,SangerSeq_start,AmplifiedSeq_start);


%keyboard
collect_number750 = [];
collect_number450 = [];
for i=1:length(inds)
  
  ampData = [];
  
  for j=1:length(inds{i})  
    
    % amplified name and number of errors - with minor
    ampData{j} = num2str([res{inds{i}(j)}{1}.Full.number length(find(align_mat(i,:)==0)) length(find(align_vec_mat(i,:)==0))]);
    ampData{j} = delSpaces(ampData{j});
    ampData{j} = ['{\color{magenta} ',ampData{j},'}'];
    data_ampData{j} = [res{inds{i}(j)}{1}.Full.number length(find(align_mat(i,:)==0)) length(find(align_vec_mat(i,:)==0))];
    
    % 750
    tmp_lab750{j} = [];
    for k=1:length(res{inds{i}(j)})
      clear lab750 lab450
      if ~isempty(res{inds{i}(j)}{k}.s750.Sequence)
        number750 = res{inds{i}(j)}{k}.s750.number(1);
        collect_number750 = [collect_number750,number750];
        data750 = [ overlapData{inds{i}(j)}{k}.s750.length overlapData{inds{i}(j)}{k}.s750.numErrorsMajor overlapData{inds{i}(j)}{k}.s750.numErrorsMinor]; 
        
        % if the bacteria was found - then boldface - otherwise - blue
        
        res750 = [zeros(3,3)];
        if ~isempty(find(tmp_ind_W_Solexa==number750))
          
          if ~isempty(find(sameNumbers==number750))
            %keyboard
            if res{inds{i}(j)}{1}.Full.number==number750
              lab750 = ['{\color{blue} {\bf "=" [',num2str(data750),']}}'];
            else
              lab750 = ['{\color{blue} {\bf ',num2str(number750),' [',num2str(data750),']}}'];
            end
            res750(1,:) = data750;
          else % only Solexa
            if res{inds{i}(j)}{1}.Full.number==number750
              lab750 = ['{\color{green} {\bf "=" [',num2str(data750),']}}'];
            else
              lab750 = ['{\color{green} {\bf ',num2str(number750),' [',num2str(data750),']}}'];
            end
            res750(2,:) = data750;
          end
        else % was not found
          if res{inds{i}(j)}{1}.Full.number==number750
            lab750 = ['{\color{orange} "=" [',num2str(data750),']}'];
          else
            lab750 = ['{\color{orange} ',num2str(number750),' [',num2str(data750),']}'];
          end
          
          res750(3,:) = data750;
        end
        
        [res{inds{i}(j)}{1}.Full.number number750]
        
      end
      
      if exist('lab750')
        tmp_lab750{j}{k} = delSpaces([lab750]);
        data_tmp_lab750{j}{k} = res750;
      else
        tmp_lab750{j}{k} = '{\color{orange} []}';
        data_tmp_lab750{j}{k} = [];
      end
      
      
    end % k
    
    
    % 450
    tmp_lab450{j} = [];
    for k=1:length(res{inds{i}(j)})
      
      if ~isempty(res{inds{i}(j)}{k}.s450.Sequence)
        %keyboard
        number450 = res{inds{i}(j)}{k}.s450.number(1);
        collect_number450 = [collect_number450,number450];
        data450 = [ overlapData{inds{i}(j)}{k}.s450.length overlapData{inds{i}(j)}{k}.s450.numErrorsMajor overlapData{inds{i}(j)}{k}.s450.numErrorsMinor]; 
        
        
        % if the bacteria was found - then boldface - otherwise - blue
        res450 = zeros(3,3);
        if ~isempty(find(tmp_ind_W_454==number450))
          if ~isempty(find(sameNumbers==number450))
            if res{inds{i}(j)}{1}.Full.number==number450
              lab450 = ['{\color{blue} {\bf "=" [',num2str(data450),']}}'];
            else
              lab450 = ['{\color{blue} {\bf ',num2str(number450),' [',num2str(data450),']}}'];
            end
            res450(1,:) = data450;
          else % only 454
            if res{inds{i}(j)}{1}.Full.number==number450
              lab450 = ['{\color{red} {\bf "=" [',num2str(data450),']}}'];
            else
              lab450 = ['{\color{red} {\bf ',num2str(number450),' [',num2str(data450),']}}'];
            end
            res450(2,:) = data450;
          end
        else % was not found
          if res{inds{i}(j)}{1}.Full.number==number450
            lab450 = ['{\color{gray} "=" [',num2str(data450),']}'];
          else
            lab450 = ['{\color{gray} ',num2str(number450),' [',num2str(data450),']}'];
          end
          res450(3,:) = data450;
        end
        
        [res{inds{i}(j)}{1}.Full.number number450]
      end
      
      if exist('lab450')
        tmp_lab450{j}{k} = delSpaces([lab450]);
        data_tmp_lab450{j}{k} = res450;
      else
        tmp_lab450{j}{k} = '{\color{gray} []}';
        data_tmp_lab450{j}{k} = [];
      end
      
      
    end
    
    
    
  end % j
  
  
  lab{i} = [];
  for j=1:length(inds{i})  
    resData{i}.ampData{j} = data_ampData{j};
    lab{i} = [lab{i},';',ampData{j}];
    for k=1:length(res{inds{i}(j)})
      resData{i}.s750{j}{k} = data_tmp_lab750{j}{k};
      resData{i}.s450{j}{k} = data_tmp_lab450{j}{k};
      lab{i} = [lab{i},';',tmp_lab750{j}{k},';',tmp_lab450{j}{k}];
    end
    
    
    %keyboard
    
  end
  
  %lab{i} = [lab{i},';'];
  %for j=1:length(inds{i})
  %  lab{i} = [lab{i},num2str(SangerOffset{i}(j)),' '];
  %end
  
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
    text(textErrors{i}.position(j),i,textErrors{i}.valAtError(:,j),'color','r','fontsize',5);
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

set(gcf,'paperorientation','landscape','paperposition',[-1    0    12    8.2]); %[from left from_bottom XX height ]
set(gca,'fontweight','bold','fontsize',fs,'position',paperPosition)

%keyboard