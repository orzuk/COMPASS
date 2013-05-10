clear
load ~/CS/BAC/12Samples/validation/scf_files_drosophila/insilocopcr_sanger_drosophila

basicDirName = '~/CS/BAC/12Samples/validation/scf_files_drosophila/';

% W forward
ampSequenceName = Sequence_amp_wolb2;
ampHeader = Header_amp_wolb2;
 = '5103_W-R.scf';
%s1_start = 60; % start of amplified sequence
%s2_start = 50; % start of Sanger major allele
%n = 600-s2_start;

%fileName = '5103_W--R.scf';
%s1_start = 60; % start of amplified sequence
%s2_start = 25; % start of Sanger major allele
%n = 720-s2_start;

fileName = '5103_W-F.scf';
s1_start = 60; % start of amplified sequence
s2_start = 30; % start of Sanger major allele
n = 620-s2_start;


if ~isempty(findstr(fileName,'F.scf'))
  F_flag = 1;
else
  F_flag = 2; % reverse
end


[seq,base_opt] = find_possible_seq_from_chromatogram_noam([basicDirName,fileName]);
[align_vec_mat,dist,xi] = compare_sangerseq_to_ampseq(ampSequenceName,n,seq,F_flag,base_opt,s1_start,s2_start,ampHeader);
 
lab = results_chromatogram_HeaderNames(ampHeader(xi),50);

% add the number of errors
for i=1:length(xi)
  lab{i} = [lab{i},' ',num2str(length(find(align_vec_mat(i,:)<1)))];
  lab{i} = strrep(lab{i},'_','\_');
end



clear tmp_ind_lab;
for i=1:length(lab)
  a = findstr(lab{i},' ');
  tmp_ind_lab(i) = str2num(lab{i}(1:a(1)-1));
end


list_Solexa_S1_S4 = load('~/CS/BAC/results_for_figure_s1_s4/listOfBAC_Solexa_S1_S4','zz');
list_454_S1_S4 = load('~/CS/BAC/results_for_figure_s1_s4/listOfBAC_454_S1_S4','zz_454');

tmp_ind_W_454 = [275390 273974 258983 52332 117492 130782 244093 91297 91789 170818];
tmp_ind_W_Solexa = [108447 564629]

sameNumbers = [108447 52332;]; % Solexa,454


% find those which appear in the amplified list and those that were not amplified

[junk,i1_Solexa,i2_Solexa] = intersect(tmp_ind_W_Solexa,tmp_ind_lab);
not_amplified_Solexa = setdiff(tmp_ind_W_Solexa,tmp_ind_lab);

[junk,i1_454,i2_454] = intersect(tmp_ind_W_454,tmp_ind_lab);
not_amplified_454 = setdiff(tmp_ind_W_454,tmp_ind_lab);

% add those not amplified
disp(not_amplified_454)



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


clf
imagesc(align_vec_mat);
colormap('gray')    

fs = 10;

%title(tit,,'fontweight','bold','fontsize',fs)

pos_xl = [100,500];
xl = cell(1,length(pos_xl))
for i=1:length(xl)
  xl{i} = num2str(pos_xl(i));
end
format_ticks(gca,xl,yl,[pos_xl],[1:size(align_vec_mat,1)])

set(gca,'fontweight','bold','fontsize',fs,'position',[0.35    0.01    0.62    0.95])

name = fileName;
name = strrep(name,'_','\_')
title(name,'fontweight','bold','fontsize',fs+20)
