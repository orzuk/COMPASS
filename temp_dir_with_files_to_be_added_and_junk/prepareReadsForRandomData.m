function [fracRelevantReadsForInfinity,normalizedBac,uniqueReads,uniqueReads_length]=prepareReadsForRandomData(readLength,Sequence1,ind_bac_in_mix,correctWeight)
disp('assumes infinite number of reads. prepareDataForRandomData')


% Sequence1 are the sequences of the current group
% load Sequence1 based on tmpInd

% seq is the list of possible sequences for each Sequence1
o = 1:length(Sequence1);
for i=1:length(o)
  if mod(i,1000)==0
    i
  end
  seq{i} = char(55*ones(length(Sequence1{o(i)})-readLength,readLength));
  for j=0:length(Sequence1{o(i)})-readLength
    seq{i}(j+1,:) = Sequence1{o(i)}(j+1:j+readLength);
  end
  seq{i} = uint8(seq{i});
end  
%keyboard
clear Sequence1 Sequence_mix

posO = [];
seqO = [];
ln = zeros(1,length(o));
for i=1:length(o)
  ln(i) = size(seq{i},1);
end
len = sum(ln);

% SeqO is a list of all sequences in all Sequence1
k = 0;
%seqO = uint8(char(55*ones(len,readLength,'uint8')));
seqO = zeros(len,readLength,'uint8');
posO = zeros(len,1,'uint16');
for i=1:length(o)
  seqO(k+1:k+length(seq{i}),:) = seq{i};
  posO(k+1:k+length(seq{i})) = i*ones(length(seq{i}),1);
  k = k+length(seq{i});
end
clear seq

%%%%%%%
% find the reads in each of these
%keyboard
% output the data
disp('assumes that the number of bacteria in part is smaller than 2^16')
%keyboard
[values,currInds,positions,leng] = myUniqueUINT8(seqO,posO);



% values = the unique rows in Esq
% currents = the indices of the unique rows in seqO
% positions = the BAC which holds the sequence and the number of times it does.
% leng = the total number of BAC which hold the sequence

clear seqO currInds posO

%keyboard
[vals_leng inds_leng] = get_duplicates(leng);
% inds_leng = indices of unique values of leng - the interesting one is leng==1

lenb = sum(leng);
clear leng

% inds_leng = indices of 

%pause
keep = 1:size(values,1);


lenKeep = length(keep);


% change to uint* ? based on number of sequences: size(values)
put = zeros(lenb,3);
k = 0;
for i=1:length(inds_leng)
  
  i1 = 1:length(inds_leng{i});

  b = inds_leng{i}(i1)*ones(1,vals_leng(i));
  
  dat = double(cell2mat(positions(inds_leng{i})));
  curr_pos = dat(1:2:end,:);curr_pos = curr_pos(i1,:);
  curr_num = dat(2:2:end,:);curr_num = curr_num(i1,:);  
  
  % curr_dat: [sequence i; bacteria j;number of times]
  curr_dat = [b(:) curr_pos(:) curr_num(:)];
  
  put(k+1:k+size(curr_dat,1),:) = curr_dat;
  k = k+size(curr_dat,1);
  
end

clear positions inds_leng d
% change this - to include values themselves and the work on the recuded number of rows


part = 1:10^5:lenb;
if part(end)<lenb
  part(end+1) = lenb;
end
part = [part,lenb+1];

normalizedBac = spalloc(length(keep),length(o),lenb);
for i=1:length(part)-1
  
  inWhichBAC = spalloc(length(keep),length(o),lenb);
  do = part(i):part(i+1)-1;
  inWhichBAC((put(do,2)-1)*lenKeep+put(do,1)) = put(do,3);
  normalizedBac = normalizedBac+inWhichBAC;
end
clear put
nnz(normalizedBac)
clear inWhichBAC keep

for i=1:length(o)
  normalizedBac(:,i) = normalizedBac(:,i)/ln(i);
end
%keyboard


% check the rank
%keyboard
uniqueReads = values;
%keyboard
%keyboard
A = spalloc(size(normalizedBac,1),size(normalizedBac,2),nnz(normalizedBac));
A(find(normalizedBac)) =1;
uniqueReads_length = sum(A,2);

if exist('ind_bac_in_mix')
  fracRelevantReadsForInfinity = normalizedBac*correctWeight';
else
  fracRelevantReadsForInfinity = 1;
end


