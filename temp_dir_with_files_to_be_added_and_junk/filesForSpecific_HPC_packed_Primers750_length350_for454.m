%pack the sequences in 64 bit
clear
load ~/CS/BAC/primers_startAs750_length350/bac16s_1to350_primers450

n = length(seq_uni1to350_primers450);
Sequence_packed64 = cell(n,1);
for i=1:n
  [Sequence_packed64{i}, len_uni] = pack_seqs(seq_uni1to350_primers450{i},64);
end

len_uni = zeros(1,n);
for i=1:n
  len_uni(i) = length(seq_uni1to350_primers450{i});
end

save('~/CS/BAC/s16_data_primers_startAs750_length350_packed64','Sequence_packed64','len_uni')    




%%%%%%%%%%%%%55
% find those which have non ACGT
clear

load ~/CS/BAC/primers_startAs750_length350/bac16s_1to350_primers450
d = [];
for i=1:length(seq_uni1to350_primers450)
  if ~isempty(find(seq_uni1to350_primers450{i}~='A' & seq_uni1to350_primers450{i}~='C'  & seq_uni1to350_primers450{i}~='G'  & seq_uni1to350_primers450{i}~='T' ))
    d = [d,i];
  end
end
save ~/CS/BAC/listOfDoubles_primers_startAs750_length350 d



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% create a new database that does not contain any non-ACGT
clear
load ~/CS/BAC/listOfDoubles_primers_startAs750_length350
load ~/CS/BAC/s16_data_primers_startAs750_length350_packed64
load('~/CS/BAC/primers_startAs750_length350/bac16s_1to350_primers450','header_uni1to350_primers450')
Sequence_packed64(d) = [];
len_uni(d) = [];
header_uni1to350_primers450(d) = [];


n = size(Sequence_packed64,1);
part = 1:100:n;
part(end) = n+1;

positionInPart = zeros(1,n);
for i=1:length(part)-1
  partition{i} = part(i):part(i+1)-1;
  positionInPart(partition{i}) = i;
end
save ~/CS/BAC/primers_startAs750_length350/datNoNonACGT/keyNoNonACGT_primers_startAs750_length350 partition part positionInPart  len_uni

for i=1:length(part)-1
  clear Seq_tmp Head_tmp 
  for j=partition{i}
    s = sprintf('seq_%d', j);
    h = sprintf('head_%d',j);
    Seq_tmp.(s) = Sequence_packed64(j);
    Head_tmp.(h) = header_uni1to350_primers450(j);
 
  end
  
  list_seq =  ['save ~/CS/BAC/primers_startAs750_length350/datNoNonACGT/packed64/seq_part_',num2str(i),' -struct Seq_tmp'];
  eval(list_seq)
  
  list_head = ['save ~/CS/BAC/primers_startAs750_length350/datNoNonACGT/packed64/head_part_',num2str(i),' -struct Head_tmp'];
  eval(list_head)
  
end

% create one dataset of the same data
clear
load ~/CS/BAC/listOfDoubles_primers_startAs750_length350
load('~/CS/BAC/primers_startAs750_length350/bac16s_1to350_primers450','header_uni1to350_primers450')

header_uni1to350_primers450(d) = [];

save ~/CS/BAC/primers_startAs750_length350/bac16s_1to350_primers450_HeaderForRunList header_uni1to350_primers450

% create one dataset of the same data - sequences and header
clear
load ~/CS/BAC/listOfDoubles_primers_startAs750_length350
load('~/CS/BAC/primers_startAs750_length350/bac16s_1to350_primers450')

header_uni1to350_primers450(d) = [];
seq_uni1to350_primers450(d) = [];
save ~/CS/BAC/primers_startAs750_length350/bac16s_1to350_primers450_SequenceAndHeader seq_uni1to350_primers450 header_uni1to350_primers450 

