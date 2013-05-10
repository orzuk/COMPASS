function [strinfilename,red]=fun_prep_data_for_simulations_loadSpecificBAC2OrFourth(basicSeqKey,basicSeqNameDir,indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag,seed_in,outfilename,N,createReadsFlag,saveToFileFlag,addNoiseFlag,ErrorStruct,ind_bac_in_mix_input)


%disp('does not save a new file name if a seed is entered. fun_prep_data_for_simulations_loadSpecificBAC2.m')
%keyboard
disp('fun_prep_data_for_simulations_loadSpecificBAC2Or.m. uses packed format')

if ~exist('addNoiseFlag','var')  
  addNoiseFlag = 0;
end
 
if addNoiseFlag 
  if ~exist('ErrorStruct','var')
    disp('ErrorStruct does not exist. fun_prep_data_for_simulations_loadSpecificBAC2Or.m');
    return
  end
end

if ~exist('saveToFileFlag','var')  
  saveToFileFlag = 0;
end

if isempty(createReadsFlag)
  createReadsFlag = 0;
end
  
%keyboard

% Sequence_uni are the database variables
% indsimbac_vec - vector of indices of similar bacteria pre-calculated for
% the case of bac_dist_flag =1.
% Nbac_in_mixture - number of bacteria in mix
% readlen - length of simulated reads in bp
% npower - power of the frequencies profile (should be greater than zero)
% bac_dist_flag - flag for the distance of the bacteria in the mix. 0 =
% random sampleing, 1 = close bacteria.
% seed_in - option to input the seed, enter [] for setting the clock as the
% seed.
% outfilename - string for file name

if isempty(seed_in)
  outseed = sum(100*clock);
else
  outseed=seed_in;  
end
rand('seed',outseed);
%keyboard


if bac_dist_flag==0
  disp('random selection')  
  ind_bac_in_mix=randi([1,N],Nbac_in_mixture,1);
else
  disp('close bacteria')
  tmp=randi([1+Nbac_in_mixture,length(indsimbac_vec)-Nbac_in_mixture],1,1);
  ind_bac_in_mix=indsimbac_vec(tmp-ceil(Nbac_in_mixture/2):tmp+ceil(Nbac_in_mixture/2)-1);
end

ind_bac_in_mix=ind_bac_in_mix(randperm(Nbac_in_mixture));


if exist('ind_bac_in_mix_input','var') && ~isempty(ind_bac_in_mix_input)
  ind_bac_in_mix = ind_bac_in_mix_input.nonZeroBAC;
  allBAC  = ind_bac_in_mix_input.allBAC;
  disp('entered ind_bac_in_mix - changed this')
  return
end


w=1./([1:Nbac_in_mixture].^npower);%setting the frequencies as the power law

if ~isinf(Nreads)
  w=round(w/sum(w)*Nreads);%this is how many read  we will have from each bacteria in the mix
                           %keyboard
  correctWeight = zeros(N,1);
  %11,correctWeight(251001:252000) = 0;
  correctWeight(ind_bac_in_mix ) = w./sum(w);
  
else
  
  w=w/sum(w);%this is how many read  we will have from each bacteria in the mix
                           %keyboard
  correctWeight = zeros(N,1);
  %11,correctWeight(251001:252000) = 0;
  correctWeight(ind_bac_in_mix ) = w;
  w = -1;
end


npowerName = num2str(npower);
npowerName(find(npowerName=='.')) = [];

basicStrintFileName = [outfilename,'_Nbacmix_',num2str(Nbac_in_mixture),'_Nread_',num2str(Nreads)...
    ,'_Readlen_',num2str(readlen),'_Npower_',npowerName,'_bacdistflag_',num2str(bac_dist_flag)];

%strinfilename=[basicStrintFileName,'_reads.mat'];
strinfilename = [basicStrintFileName,'_NoReads.mat'];


%keyboard

if isempty(seed_in) | saveToFileFlag==1 % if no seed - saves to file
  save(strinfilename,'Nbac_in_mixture', 'ind_bac_in_mix', 'w', 'outseed','correctWeight','Nreads','readlen','npower','bac_dist_flag');
end



if createReadsFlag==0
  red = [];
  return
end
clear correctWeight



load(basicSeqKey,'positionInPart','len_uni');
len_mix = len_uni(ind_bac_in_mix); clear len_uni
Sequence_mix = cell(length(ind_bac_in_mix),1);
for i=1:length(ind_bac_in_mix)
  clear seq_*
  load([basicSeqNameDir,'seq_part_',num2str(positionInPart(ind_bac_in_mix(i)))],['seq_',num2str(ind_bac_in_mix(i))]);
  ww = ['Sequence_mix{i} = seq_',num2str(ind_bac_in_mix(i)),'{1};'];
  eval(ww);
  
end
clear seq_* positionInPart head_*


%keyboard

%%%%%%%%%5
% old version
if 1==2
  
  
  %rand('seed',1)  
  bac_read_startvec=zeros(Nreads,1);
  k=1;
  for i=1:Nbac_in_mixture
    bac_read_startvec(k:k+w(i)-1)=randi([1,len_mix(i)-readlen+1],w(i),1);
    k=k+w(i);
  end
  
  currSeq = extract_sub_kmers(Sequence_mix{1}, len_mix(1), readlen, 0,0,1,bac_read_startvec(1));
  red = zeros(Nreads,length(currSeq),'uint64');
  % end change
  
  %maxLength = 0;
 
  k=1;
  % change Nov 2011
  for i=1:Nbac_in_mixture
    for j=1:w(i)
      currSeq = extract_sub_kmers(Sequence_mix{i}, len_mix(i), readlen, 0,0,1,bac_read_startvec(k+j-1));
      %currN = length(currSeq);
      %maxLength = max([maxLength,currN]);                           
      %red(k+j-1,1:currN) = currSeq;
      red(k+j-1,:) = currSeq;
    end
    k=k+w(i);
  end
  clear bac_read_startvec

  red(k:end,:) = [];
  
end
%%%%%%%%%%%%%%%%%%%%%%%% end old version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%keyboard

currSeq = extract_sub_kmers(Sequence_mix{1}, len_mix(1), readlen, 0,0,1,1);
red = zeros(Nreads,length(currSeq),'uint64');clear currSeq
k = 1;
for i=1:Nbac_in_mixture
  bac_read_startvec=randi([1,len_mix(i)-readlen+1],w(i),1);
  for j=1:w(i)
    red(k+j-1,:) = extract_sub_kmers(Sequence_mix{i}, len_mix(i), readlen, 0,0,1,bac_read_startvec(j));
  end
  
  k=k+w(i);  
end

%return

%red(:,maxLength+1:end) = [];

% add noise to red

if addNoiseFlag
  read_error_substitution_table = GenerateSubstitutionErrorTable(readlen, ErrorStruct);
  
  %red = add_noise_to_kmers(readlen, size(red,1), ...
  %      red, reshape(read_error_substitution_table, readlen, 4*4));
  if size(red,1)>10^5
    partRed = 1:10^5:size(red,1);
    partRed(end+1) = size(red,1)+1;
    
    for ind_part_red=1:length(partRed)-1
      curPartRed = partRed(ind_part_red):partRed(ind_part_red+1)-1;
      red(curPartRed,:) = add_noise_to_kmers(readlen, size(curPartRed,2), ...
        red(curPartRed,:), reshape(read_error_substitution_table, readlen, 4*4));
    end
  end
  
  %keyboard
  disp(['added noised according to: model: ',ErrorStruct.error_model,' baseline_error: ',num2str(ErrorStruct.baseline_error),' final_error: ',num2str(ErrorStruct.final_error)]);

end
%keyboard
%compare the former and the current - compared: see fun_prep_data_for_simulations_loadSpecificBAC_checkOrAndFormer

% what does he do when length are different

tmp=randperm(size(red,1));
red=red(tmp,:);
%keyboard