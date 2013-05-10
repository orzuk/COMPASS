%split 750 450  - what is the order we make? open names

function [errorsMinorOverSequence,reordered_cData]=results_chroma6_in(basicDirName,fileName,cData,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,tmp_ind_W_454,tmp_ind_W_Solexa,sameNumbers,paperPosition,thresh,fs,th_sec2first)
%keyboar
lenName = 50;

if ~isempty(findstr(fileName,'F.scf'))
  F_flag = 1;
else
  F_flag = 2; % reverse
end


% find Sanger sequences
[seqSanger,base_opt] = find_possible_seq_from_chromatogram_noam3([basicDirName,fileName],th_sec2first);

% collect data and errors and overlap. note it's reordered,. add offset add errors over sequence
[align_mat,align_vec_mat,overlapData,saveErrors,SangerOffset,reordered_cData,errorsMinorOverSequence]=compare_sangerseq_to_ampseq4(cData,lengthValidSanger,seqSanger,F_flag,base_opt,SangerSeq_start,AmplifiedSeq_start);

% find those that were found - all their names
found = struct;
%keyboard
%found.s750 = [tmp_ind_W_Solexa',tmp_ind_W_Solexa'];
found.s750 = [];
for i=1:length(cData.s750)
  flag = 0;
  currNum = zeros(length(cData.s750{i}.Header),1);
  for j=1:length(cData.s750{i}.Header)
    currNum(j) = findNumber(cData.s750{i}.Header{j});
  end
  
  int = intersect(currNum,tmp_ind_W_Solexa);
  if ~isempty(int)
    k = size(found.s750,1);
    found.s750(k+1:k+length(currNum),1) = [currNum];
    found.s750(k+1:k+length(currNum),2) = int;
  end
end

%found.s450 = [tmp_ind_W_454',tmp_ind_W_454'];
found.s450 = [];
for i=1:length(cData.s450)
  flag = 0;
  currNum = zeros(length(cData.s450{i}.Header),1);
  for j=1:length(cData.s450{i}.Header)
    currNum(j) = findNumber(cData.s450{i}.Header{j});
  end
  %if ~isempty(find(currNum==275390)),i,pause;end
  int = intersect(currNum,tmp_ind_W_454);
  if ~isempty(int)
    k = size(found.s450,1);
    found.s450(k+1:k+length(currNum),1) = [currNum];
    found.s450(k+1:k+length(currNum),2) = int;
  end
end


[junk,i1] =  unique(found.s750(:,1));
found.s750 = found.s750(i1,:);
[junk,i1] =  unique(found.s450(:,1));
found.s450 = found.s450(i1,:);


% find those that were not amplified and appear in our list
not_amplified_Solexa = setdiff(tmp_ind_W_Solexa,found.s750(:,1));
not_amplified_454 = setdiff(tmp_ind_W_454,found.s450(:,1));
%keyboard

%%%%%%%%
% plot
clear lab
for i=1:length(reordered_cData.Numbers)
  currNumber = reordered_cData.Numbers(i);
  lab{i}{1} = sprintf('%s\t %s\t %s\t %s',num2str(currNumber),num2str(overlapData.amp(i,1)),num2str(overlapData.amp(i,2)),num2str(overlapData.amp(i,3)));
  %['{\color{magenta}[',num2str(currNumber),' ',num2str(overlapData.amp(i,:)),']};'];
  [junk,ind] = intersect(found.s750(:,1),currNumber);
  if ~isempty(ind)
    lab{i}{2} = sprintf('b_%s\t %s\t %s',num2str(found.s750(ind,2)),num2str(overlapData.s750(i,[1 3])));
              %'{\colorbox{blue}{\color{blue}}{[',num2str(found.s750(ind,2)),' ',num2str(overlapData.s750(i,:)),']}};'];
  else
    lab{i}{2} = sprintf('g_%s\t %s\t %s',...
                        num2str(overlapData.s750(i,1)),...
                        num2str(overlapData.s750(i,2)),...
                        num2str(overlapData.s750(i,3)) );
    
    %[lab{i},'{\color{gray} [',num2str(overlapData.s750(i,:)),']};'];
  end
  [junk,ind] = intersect(found.s450(:,1),currNumber);
  if ~isempty(ind)
    lab{i}{3} = sprintf('g_%s\t %s\t %s\t %s',...
                        num2str(found.s450(ind,2)),...
                        num2str(overlapData.s450(i,1)),...
                        num2str(overlapData.s450(i,2)),...
                        num2str(overlapData.s450(i,3)) );
                       
    %lab{i} = [lab{i},'{\colorbox{red}{\color{red}}{[',num2str(found.s450(ind,2)),' ',num2str(overlapData.s450(i,:)),']}};'];
  else
    lab{i}{3} = sprintf('g_%s\t %s\t %s',...
                        num2str(overlapData.s450(i,1)),...
                        num2str(overlapData.s450(i,2)),...
                        num2str(overlapData.s450(i,3)) );
    %lab{i} = [lab{i},'{\color{gray} [',num2str(overlapData.s450(i,:)),']};'];
  end
  
  %lab{i} = delSpaces(lab{i});
end
% add those that were not amplified

%keyboard
% text errors
clear textErrors
for i=1:length(reordered_cData.Numbers)
  findErrors = find(saveErrors{i}(1,:)~='P');
  textErrors{i}.position = findErrors;
  for k=1:length(findErrors)
      currPos = findErrors(k)+SangerOffset(i)-1;
      if currPos>=100
        positionInText = num2str(currPos)';
      elseif currPos>=10 & currPos<100
        positionInText = [0;num2str(currPos)'];
      else % <10
        positionInText = [0;0;num2str(currPos)'];
      end
      textErrors{i}.valAtError(:,k) = [saveErrors{i}(:,findErrors(k));positionInText];
    end
end



align_vec_mat = [align_vec_mat;2*ones(length(not_amplified_Solexa),size(align_vec_mat,2))];
add_yl_Solexa = cell(1,length(not_amplified_Solexa));
for i=1:length(not_amplified_Solexa)
  add_yl_Solexa{i} = sprintf('b_%s\t NA',num2str(not_amplified_Solexa(i)));
  %add_yl_Solexa{i} = ['{\color{blue} ',num2str(not_amplified_Solexa(i)),' NA}'];
end

align_vec_mat = [align_vec_mat;2*ones(length(not_amplified_454),size(align_vec_mat,2))];
add_yl_454 = cell(1,length(not_amplified_454));
for i=1:length(not_amplified_454)
  add_yl_454{i} = sprintf('b_%s\t NA',num2str(not_amplified_454(i)));
  %add_yl_454{i} = ['{\color{red} ',num2str(not_amplified_454(i)),' NA}'];
end

xl = cell(1,length(pos_xl))
for i=1:length(xl)
  xl{i} = num2str(pos_xl(i));
end
%keyboard

keyboard
clf
subplot(1,2,2)
imagesc(align_vec_mat);
for i=1:size(align_vec_mat,1)
  line([1 size(align_vec_mat,2)],[i+0.5,i+0.5]','color','k');        
end
colormap('gray')    


%for i=1:length(textErrors)
%  for j=1:length(textErrors{i}.position)
%    text(textErrors{i}.position(j),i,textErrors{i}.valAtError(:,j),'color','r','fontsize',8);
%  end
%end



toPrint = [lab,add_yl_454,add_yl_Solexa];%



subplot(1,2,1)
printChromatogramResults(toPrint)

%if ~exist('fs')
%  fs = 10;
%end
%format_ticks(gca,xl,toPrint,[pos_xl],[1:size(align_vec_mat,1)],[],[],[],'fontweight','bold','fontsize',fs)


name = fileName;
name = strrep(name,'_','\_')
title(name,'fontweight','bold','fontsize',fs+10)
%keyboard

set(gcf,'paperorientation','landscape','paperposition',[-1    0    12    8.5]); %[from left from_bottom XX height ]
set(gca,'fontweight','bold','fontsize',fs,'position',paperPosition)

%keyboard


function currNumber=findNumber(Header)

currNumber = Header;
a = find(currNumber==' ');
currNumber = currNumber(1:a(1)-1);
currNumber = str2num(currNumber);
