%pack the sequences in 64 bit
clear
load ~/CS/BAC/bacteria_s16_data_uni

n = length(Sequence_uni);
Sequence_packed64 = cell(n,1);
for i=1:n
  [Sequence_packed64{i}, seqs_len] = pack_seqs(Sequence_uni{i},64);
end

save('~/CS/BAC/s16_data_uni_packed64','Sequence_packed64','len_uni')    


% lengths

%%%%%%%%%%%%%%%%%%%% for database - done only once
clear
load ~/CS/BAC/s16_data_uni_packed64

n = size(Header_uni,1);
part = 1:100:n;
part(end) = n+1;

positionInPart = zeros(1,n);
for i=1:length(part)-1
  partition{i} = part(i):part(i+1)-1;
  positionInPart(partition{i}) = i;
end
%save ~/CS/BAC/dat450000/key450000 partition part positionInPart  len_uni 

for i=1:length(part)-1
  clear Seq_tmp Head_tmp 
  for j=partition{i}
    s = sprintf('seq_%d', j);
    h = sprintf('head_%d',j);
    Seq_tmp.(s) = Sequence_packed64(j);
    Head_tmp.(h) = Header_uni(j);
    
  end
  list_seq = ['save ~/CS/BAC/dat450000/packed64/seq_part_',num2str(i),' -struct Seq_tmp'];
  eval(list_seq)
  list_head = ['save ~/CS/BAC/dat450000/packed64/head_part_',num2str(i),' -struct Head_tmp'];
  eval(list_head)
  
end
%%%%%%%%%%%%%%%%%% end for database - done only once


%profile -memory on;
%setpref('profiler','showJitLines',1);
%[normalizedBac values] = prepareGroupOf1000DistributedSequenceFiles(readLength,tmpInd{i},basicSeqNameDir,basicSeqKey);
%profile report;profile off

%%%%%%%%%%%%%55
% find those which have non ACGT
load ~/CS/BAC/bacteria_s16_data_uni

d = [];
for i=1:length(Sequence_uni)
  if ~isempty(find(Sequence_uni{i}~='A' & Sequence_uni{i}~='C'  & Sequence_uni{i}~='G'  & Sequence_uni{i}~='T' ))
    d = [d,i];
  end
end
save ~/CS/BAC/listOfDoubles d



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% create a new database that does not contain any non-ACGT
clear
load ~/CS/BAC/listOfDoubles
load ~/CS/BAC/s16_data_uni_packed64
load('~/CS/BAC/bacteria_s16_data_uni','Header_uni')
Sequence_packed64(d) = [];
len_uni(d) = [];
Header_uni(d) = [];


n = size(Sequence_packed64,1);
part = 1:100:n;
part(end) = n+1;

positionInPart = zeros(1,n);
for i=1:length(part)-1
  partition{i} = part(i):part(i+1)-1;
  positionInPart(partition{i}) = i;
end
save ~/CS/BAC/datNoNonACGT/keyNoNonACGT partition part positionInPart  len_uni 

for i=1:length(part)-1
  clear Seq_tmp Head_tmp 
  for j=partition{i}
    s = sprintf('seq_%d', j);
    h = sprintf('head_%d',j);
    Seq_tmp.(s) = Sequence_packed64(j);
    Head_tmp.(h) = Header_uni(j);
 
  end
  %list_seq = ['save ~/CS/BAC/datNoNonACGT/packed64/seq_part_',num2str(i),' -struct Seq_tmp'];
  %eval(list_seq)

  list_head = ['save ~/CS/BAC/datNoNonACGT/packed64/head_part_',num2str(i),' -struct Head_tmp'];
  eval(list_head)

end

% do the same for the unpacked
clear
load ~/CS/BAC/listOfDoubles
load ~/CS/BAC/bacteria_s16_data_uni

Sequence_uni(d) = [];
Header_uni(d) = [];


load ~/CS/BAC/datNoNonACGT/keyNoNonACGT
for i=1:length(part)-1
  clear Seq_tmp Head_tmp 
  for j=partition{i}
    s = sprintf('seq_%d', j);
    h = sprintf('head_%d',j);
    Seq_tmp.(s) = Sequence_uni(j);
    Head_tmp.(h) = Header_uni(j);
  end
  list_seq = ['save ~/CS/BAC/datNoNonACGT/unpacked/seq_part_',num2str(i),' -struct Seq_tmp'];
  eval(list_seq)
  
end

% create one database of the same data 
% done on 28.3.12
clear
% copied listOfDoubles from ~/CS/BAC/data_before_right_edge_correction/ to ~/CS/BAC/full16S/

load ~/CS/BAC/data_before_right_edge_correction/listOfDoubles
load ~/CS/BAC/bacteria_s16_data_uni

Header_uni(d) = [];
Sequence_uni(d) = [];
save ~/CS/BAC/full16S/bac16s_full_without_ambiguous Header_uni Sequence_uni


for i=1:size(Header_uni,1)
  if ~isempty(findstr(Header_uni{i},'161334'))
    i
  end
end


seq = 'GTGACGGTAGCTTACCAGAAAGGGACGGCTAACTACGTGCNNGCAGCCGCGGTAATACGTAGNNCCCGAGCGTNGTCCNGNATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTTGATAAGTCTGAAGTTAAAGGCTGTGGCTCAACCATAGTTCGCTTTGGAAACTGTCAAACTTGAGTGCAGAAGGGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTNTCTGGTCTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAGGTGTTGGATCCTTTCCGGGATTCAGTGCCGCAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGACCGCAAGGTTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCCGATGCTATTTCTAGAGATAGAAAGTTACTTCGGTACATCGGTGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCCTATTGTTAGTTGCCATCATTCAGTTGGGCACTCTAGCGAGACTGCCGGTAATAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGACCTGGGCTACACACGTGCTACAATGGTTGGTACAACGAGTTGCGAGTCGGTGACGGCAAGCTAATCTCTTAAAGCCAATCTCAGTTCGGATTGTAGGCTGCAACTCGCCTACATGAAGTCGGAATCGCTAGTAATCGCGGATCAGCACGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAGAGTTTGTAACACCCGAAG';

seq = 'ATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTTGATAAGTCTGAAGTTAAAGGCTGTGGCTCAACCATAGTTCGCTTTGGAAACTGTCAAACTTGAGTGCAGAAGGGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGTCTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAGGTGTTGGATCCTTTCCGGGATTCAGTGCCGCAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGACCGCAAGGTTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCCGATGCTATTTCTAGAGATAGAAAGTTACTTCGGTACATCGGTGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCCTATTGTTAGTTGCCATCATTCAGTTGGGCACTCTAGCGAGACTGCCGGTAATAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGACCTGGGCTACACACGTGCTACAATGGTTGGTACAACGAGTTGCGAGTCGGTGACGGCAAGCTAATCTCTTAAAGCCAATCTCAGTTCGGATTGTAGGCTGCAACTCGCCTACATGAAGTCGGAATCGCTAGTAATCGCGGATCAGCACGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAGAGTTTGTAACACCCGAAG';


userDir = getuserdir
load([userDir,'/CS/BAC/12Samples/Solexa/data/forAlign_full16S/dataForAlign_full16S']) 

get_mean_posalign_rna(seq,database1,database_size,bacseq_len)


