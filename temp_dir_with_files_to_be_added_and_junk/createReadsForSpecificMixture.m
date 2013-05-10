function [red]=createReadsForSpecificMixture(auxData)
%keyboard
%keyboard
if ~isfield(auxData,'seed_in')
  outseed = sum(100*clock);
else
  outseed = auxData.seed_in;  
end
rand('seed',outseed);

correctWeight = auxData.correctWeight;
ind_bac_in_mix = find(correctWeight);
Nreads = auxData.Nreads;

w = correctWeight(ind_bac_in_mix);
w=round(w/sum(w)*Nreads);%


load(auxData.basicSeqKey,'positionInPart','len_uni');
len_mix = len_uni(ind_bac_in_mix); clear len_uni
Sequence_mix = cell(length(ind_bac_in_mix),1);
for i=1:length(ind_bac_in_mix)
  clear seq_*
  load([auxData.basicSeqNameDir,'seq_part_',num2str(positionInPart(ind_bac_in_mix(i)))],['seq_',num2str(ind_bac_in_mix(i))]);
  ww = ['Sequence_mix{i} = seq_',num2str(ind_bac_in_mix(i)),'{1};'];
  eval(ww);
  
  %clear head_*
  %load([auxData.basicSeqNameDir,'head_part_',num2str(positionInPart(ind_bac_in_mix(i)))],['head_',num2str(ind_bac_in_mix(i))]);
  %ww = ['Header_mix{i} = head_',num2str(ind_bac_in_mix(i)),'{1};'];
  %eval(ww);
  
end
clear seq_* positionInPart 


Nbac_in_mixture = length(w);
currSeq = extract_sub_kmers(Sequence_mix{1}, len_mix(1), auxData.readLength, 0,0,1,1);
red = zeros(Nreads,length(currSeq),'uint64');clear currSeq
k = 1;
%keyboard
for i=1:Nbac_in_mixture
  bac_read_startvec=randi([1,len_mix(i)-auxData.readLength+1],w(i),1);
  for j=1:w(i)
    red(k+j-1,:) = extract_sub_kmers(Sequence_mix{i}, len_mix(i), auxData.readLength, 0,0,1,bac_read_startvec(j));
  end
  k=k+w(i);  
end

% change 22.6.12
red(k:end,:) = [];

if auxData.addNoiseFlag
  read_error_substitution_table = GenerateSubstitutionErrorTable(auxData.readLength, auxData.ErrorStruct);
  
  %red = add_noise_to_kmers(auxData.readLength, size(red,1), ...
  %      red, reshape(read_error_substitution_table, auxData.readLength, 4*4));
  if size(red,1)>10^5
    partRed = 1:10^5:size(red,1);
    partRed(end+1) = size(red,1)+1;
    
    for ind_part_red=1:length(partRed)-1
      curPartRed = partRed(ind_part_red):partRed(ind_part_red+1)-1;
      red(curPartRed,:) = add_noise_to_kmers(auxData.readLength, size(curPartRed,2), ...
        red(curPartRed,:), reshape(read_error_substitution_table, auxData.readLength, 4*4));
    end
  else % for less than 100000 - one unit of noise
    %keyboard
    partRed = 1;
    partRed(end+1) = size(red,1)+1;
    for ind_part_red=1:length(partRed)-1
      curPartRed = partRed(ind_part_red):partRed(ind_part_red+1)-1;
      red(curPartRed,:) = add_noise_to_kmers(auxData.readLength, size(curPartRed,2), ...
        red(curPartRed,:), reshape(read_error_substitution_table, auxData.readLength, 4*4));
    end
    
  end
  
  
  %keyboard
  disp(['added noised according to: model: ',auxData.ErrorStruct.error_model,' baseline_error: ',num2str(auxData.ErrorStruct.baseline_error),' final_error: ',num2str(auxData.ErrorStruct.final_error)]);

end

tmp=randperm(size(red,1));
red=red(tmp,:);

