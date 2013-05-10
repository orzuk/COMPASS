function [strinfilename,red]=fun_prep_data_for_simulations_loadSpecificBAC2_checkRegAndPacked(basicSeqKey,basicSeqNameDir,indsimbac_vec,Nbac_in_mixture,Nreads,readlen,npower,bac_dist_flag,seed_in,outfilename,N,createReadsFlag,saveToFileFlag,ind_bac_in_mix_input)


%disp('does not save a new file name if a seed is entered. fun_prep_data_for_simulations_loadSpecificBAC2.m')
%keyboard
disp('fun_prep_data_for_simulations_loadSpecificBAC2Or.m. uses packed format')

if ~exist('saveToFileFlag')  saveToFileFlag = 0;
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
%keyboard

if exist('ind_bac_in_mix_input')
  ind_bac_in_mix = ind_bac_in_mix_input.nonZeroBAC;
  allBAC  = ind_bac_in_mix_input.allBAC;
  disp('entered ind_bac_in_mix - changed this')
  pause
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



%keyboard

if 1==1
  load(basicSeqKey,'positionInPart','len_uni')
  Sequence_mix = cell(length(ind_bac_in_mix),1);
  for i=1:length(ind_bac_in_mix)
    clear seq_*
    load([basicSeqNameDir,'seq_part_',num2str(positionInPart(ind_bac_in_mix(i)))],['seq_',num2str(ind_bac_in_mix(i))]);
    ww = ['Sequence_mix{i} = seq_',num2str(ind_bac_in_mix(i)),'{1};'];
    eval(ww);
    
  end
  clear seq_* positionInPart head_*


  len_mix = len_uni(ind_bac_in_mix);
  bac_read_startvec=zeros(Nreads,1);
  k=1;
  for i=1:Nbac_in_mixture
    bac_read_startvec(k:k+w(i)-1)=randi([1,len_mix(i)-readlen+1],w(i),1);
    k=k+w(i);
  end

  keyboard

  seq_lens_in_words = length_cell(Sequence_mix);
  red = zeros(Nreads,max(seq_lens_in_words),'uint64');
  maxLength = 0;
  k=1;
  for i=1:Nbac_in_mixture
    i  
    for j=1:w(i)
      currSeq = extract_sub_kmers(Sequence_mix{i}, len_mix(i), readlen, 0,0,1,bac_read_startvec(k+j-1));
      currN = length(currSeq);
      maxLength = max([maxLength,currN])                           
      red(k+j-1,1:currN) = currSeq;
      [i,j]
    end
    k=k+w(i);
  end

  red(k:end,:) = [];
  red(:,maxLength+1:end) = [];

else % unpacked
  load(basicSeqKey,'positionInPart','len_uni')

  Sequence_mix1 = cell(length(ind_bac_in_mix),1);
  for i=1:length(ind_bac_in_mix)
    clear seq_*
    load(['/home/csfaculty/shental/CS/BAC/dat450000/','seq_part_',num2str(positionInPart(ind_bac_in_mix(i)))],['seq_',num2str(ind_bac_in_mix(i))]);
    ww = ['Sequence_mix1{i} = seq_',num2str(ind_bac_in_mix(i)),'{1};'];
    eval(ww);
  end
  clear seq_* positionInPart head_*

  %keyboard

  red1=repmat('n',Nreads,readlen);
  k=1;
  %keyboard
  for i=1:Nbac_in_mixture
    for j=1:w(i)
      red1(k+j-1,:)=Sequence_mix1{i}(bac_read_startvec(k+j-1):bac_read_startvec(k+j-1)+readlen-1);
    end
    k=k+w(i);
  end
  red1(k:end,:) = [];
  
  
  
  for i=200000:size(red1,1)
    a = unpack_seqs(red(i,:),50,64);a = int2nt(a);
    b = find(red1(i,:)-a);
    if ~isempty(b) && ~isempty(find(red1(i,b)=='A' | red1(i,b)=='C' | red1(i,b)=='G' | red1(i,b)=='T'))
      
      i
      pause
    end
  end
    
end



%compare the former and the current
% what does he do when length are different

tmp=randperm(size(red,1));
red=red(tmp,:);
