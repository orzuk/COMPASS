% pack the sequences in 64 bit
clear
databaseDir = ['~/CS/BAC/data750withMM/'];
databaseName = ['database_16s_750primers_with4mm']
load([databaseDir,databaseName]) 

% Sequence_uni_amp_uni(1) =[]; Header_uni_amp_uni(1) = [];

n = length(Sequence_uni_amp_uni);
Sequence_packed64 = cell(n,1);
for i=1:n
  [Sequence_packed64{i}, len_uni] = pack_seqs(Sequence_uni_amp_uni{i},64);
end

len_uni = zeros(1,n);
for i=1:n
  len_uni(i) = length(Sequence_uni_amp_uni{i});
end

save([databaseDir,databaseName,'_packed64'],'Sequence_packed64','len_uni')    


%%%%%%%%%%%%%55
% find those which have non ACGT
clear
databaseDir = ['~/CS/BAC/data750withMM/'];
databaseName = ['database_16s_750primers_with4mm']


load([databaseDir,databaseName])
d = [];
for i=1:length(Sequence_uni_amp_uni)
  if ~isempty(find(Sequence_uni_amp_uni{i}~='A' & Sequence_uni_amp_uni{i}~='C'  & Sequence_uni_amp_uni{i}~='G'  & Sequence_uni_amp_uni{i}~='T' ))
    d = [d,i];
  end
end
save([databaseDir,databaseName,'_listOfDoubles'],'d') 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% create a new database that does not contain any non-ACGT
clear
databaseDir = ['~/CS/BAC/data750withMM/'];
databaseName = ['database_16s_750primers_with4mm']

load([databaseDir,databaseName,'_listOfDoubles']) 
load([databaseDir,databaseName,'_packed64']) 
load([databaseDir,databaseName],'Header_uni_amp_uni')

Sequence_packed64(d) = [];
len_uni(d) = [];
Header_uni_amp_uni(d) = [];


n = size(Sequence_packed64,1);
part = 1:100:n;
part(end) = n+1;

positionInPart = zeros(1,n);
for i=1:length(part)-1
  partition{i} = part(i):part(i+1)-1;
  positionInPart(partition{i}) = i;
end
save([databaseDir,'datNoNonACGT/keyNoNonACGT_',databaseName],'partition','part','positionInPart','len_uni');  

for i=1:length(part)-1
  clear Seq_tmp Head_tmp 
  for j=partition{i}
    s = sprintf('seq_%d', j);
    h = sprintf('head_%d',j);
    Seq_tmp.(s) = Sequence_packed64(j);
    Head_tmp.(h) = Header_uni_amp_uni(j);
 
  end
  
  list_seq =  ['save ',databaseDir,'datNoNonACGT/packed64/seq_part_',num2str(i),' -struct Seq_tmp'];
  eval(list_seq)
  
  list_head = ['save ',databaseDir,'datNoNonACGT/packed64/head_part_',num2str(i),' -struct Head_tmp'];
  eval(list_head)
  
end


% create one database of the same data 
clear
databaseDir = ['~/CS/BAC/data750withMM/'];
databaseName = ['database_16s_750primers_with4mm']

load([databaseDir,databaseName,'_listOfDoubles']) 
load([databaseDir,databaseName])

Sequence_no_nonACGT = Sequence_uni_amp_uni;
Header_no_nonACGT = Header_uni_amp_uni;

Header_no_nonACGT(d) = [];
Sequence_no_nonACGT(d) = [];



save([databaseDir,databaseName,'_full_without_ambiguous'],'Sequence_no_nonACGT','Header_no_nonACGT')



