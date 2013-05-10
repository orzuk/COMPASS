function results_chroma_in(basicDirName,fileName,ampSequenceName,ampHeader,s1_start,s2_start,n,pos_xl,tmp_ind_W_454,tmp_ind_W_Solexa,list_454_S1_S4,list_Solexa_S1_S4,sameNumbers)

lenName = 50;

if ~isempty(findstr(fileName,'F.scf'))
  F_flag = 1;
else
  F_flag = 2; % reverse
end
%keyboard
[seq,base_opt] = find_possible_seq_from_chromatogram_noam([basicDirName,fileName]);
% dd the classSeq - chnage this
[align_vec_mat,dist,xi,classSeqs] = compare_sangerseq_to_ampseq(ampSequenceName,n,seq,F_flag,base_opt,s1_start,s2_start,ampHeader);
 
lab = results_chromatogram_HeaderNames(ampHeader(xi),lenName);

% add the number of errors
for i=1:length(xi)
  lab{i} = [lab{i},' ',num2str(classSeqs(i)),' ',num2str(length(find(align_vec_mat(i,:)<1)))];
  lab{i} = strrep(lab{i},'\','\\');
  lab{i} = strrep(lab{i},'_','\_');
end


clear tmp_ind_lab;
for i=1:length(lab)
  a = findstr(lab{i},' ');
  tmp_ind_lab(i) = str2num(lab{i}(1:a(1)-1));
end

%keyboard
[junk,i1_Solexa,i2_Solexa] = intersect(tmp_ind_W_Solexa,tmp_ind_lab);
not_amplified_Solexa = setdiff(tmp_ind_W_Solexa,tmp_ind_lab);

[junk,i1_454,i2_454] = intersect(tmp_ind_W_454,tmp_ind_lab);
not_amplified_454 = setdiff(tmp_ind_W_454,tmp_ind_lab);

add_Header_not_amplified_Solexa = cell(length(not_amplified_Solexa),1);
for i=1:length(not_amplified_Solexa)
  flag = 0;
  j = 1;
  while flag==0 && j<=size(list_Solexa_S1_S4.zz,1)
    name = list_Solexa_S1_S4.zz{j,1};
    a = find(name==' ');
    name = name(1:a(1)-1);
    name = str2num(name);
    if name==not_amplified_Solexa(i)
      flag = 1;
    else
      j = j+1;
    end
  end
  if j>size(list_Solexa_S1_S4.zz,1)
    disp('problem')
  end
  add_Header_not_amplified_Solexa{i} = list_Solexa_S1_S4.zz{j,1};
end

add_Header_not_amplified_454 = cell(length(not_amplified_454),1);
for i=1:length(not_amplified_454)
  flag = 0;
  j = 1;
  while flag==0
    name = list_454_S1_S4.zz_454{j,1};
    a = find(name==' ');
    name = name(1:a(1)-1);
    name = str2num(name);
    if name==not_amplified_454(i)
      flag = 1;
    else
      j = j+1;
    end
  end
  if j>size(list_454_S1_S4.zz_454,1)
    disp('problem')
  end
  add_Header_not_amplified_454{i} = list_454_S1_S4.zz_454{j,1};
end





% find those which appear in the amplified list and those that were not amplified




%%%%%%%%
% plot

yl = cell(1,size(align_vec_mat,1));
for i=1:size(align_vec_mat,1)
  if ~isempty(find(i2_Solexa==i)) & isempty(find(i2_454==i)) 
    yl{i} = ['{\color{green} ',lab{i},'}'];
  end
  
  if isempty(find(i2_Solexa==i)) & ~isempty(find(i2_454==i))
    yl{i} = ['{\color{red} ',lab{i},'}'];
  end

  % overwrites in case they are the same
  if (~isempty(find(i2_Solexa==i)) & ~isempty(find(i2_454==i))) | ...
        (~isempty(find(sameNumbers==tmp_ind_lab(i))))
    yl{i} = ['{\color{blue} ',lab{i},'}'];
  end
  
  if isempty(find(i2_Solexa==i)) & isempty(find(i2_454==i)) % not found by both
    yl{i} = ['{\color{gray} ',lab{i},'}'];
  end
end



% add those that were not amplified

align_vec_mat = [align_vec_mat;2*ones(length(add_Header_not_amplified_454),size(align_vec_mat,2))];
%keyboard
add_yl_454 = cell(1,length(add_Header_not_amplified_454));
for i=1:length(add_Header_not_amplified_454)
  curr_lenName = length(add_Header_not_amplified_454{i});
  add_yl_454{i} = ['{\color{red} ',add_Header_not_amplified_454{i}(1:min(curr_lenName,lenName)),' NA}'];
end

align_vec_mat = [align_vec_mat;2*ones(length(add_Header_not_amplified_Solexa),size(align_vec_mat,2))];

add_yl_Solexa = cell(1,length(add_Header_not_amplified_Solexa));
for i=1:length(add_Header_not_amplified_Solexa)
  curr_lenName = length(add_Header_not_amplified_Solexa{i});
  add_yl_Solexa{i} = ['{\color{cyan} ',add_Header_not_amplified_Solexa{i}(1:min(curr_lenName,lenName)),' NA}'];
end



clf
imagesc(align_vec_mat);
colormap('gray')    


xl = cell(1,length(pos_xl))
for i=1:length(xl)
  xl{i} = num2str(pos_xl(i));
end
%keyboard

fs = 6;
format_ticks(gca,xl,[yl,add_yl_454,add_yl_Solexa],[pos_xl],[1:size(align_vec_mat,1)],[],[],[],'fontweight','bold','fontsize',fs)

set(gca,'fontweight','bold','fontsize',fs,'position',[0.35    0.01    0.62    0.91])

name = fileName;
name = strrep(name,'_','\_')
title(name,'fontweight','bold','fontsize',fs+10)

