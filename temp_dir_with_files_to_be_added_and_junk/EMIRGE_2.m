clear

numBac = 410849;
bact = ceil(rand(1,200)*numBac);
auxData = struct;
auxData.correctWeight = zeros(1,numBac);
auxData.correctWeight(bact) = 1/length(bact);
auxData.Nreads = 10^6;
auxData.basicSeqKey = '~/CS/BAC/datNoNonACGT/keyNoNonACGT';
auxData.basicSeqNameDir = '~/CS/BAC/datNoNonACGT/packed64/';
auxData.addNoiseFlag = 0;
auxData.readLength = 100;
[red]=createReadsForSpecificMixture(auxData);

red = int2nt(unpack_seqs(red,auxData.readLength*ones(size(red,1),1),64));

load('~/CS/BAC/full16S/bac16s_full_without_ambiguous')

fileName = '200';
save(['~/CS/BAC/EMIRGE/csmiller-EMIRGE-85bf31a/try/listOfBact_',fileName],'bact');


% create reads
reg_1 = red;

clear data_1
data_1 = struct;
for i=1:size(reg_1,1)
  data_1(i) = struct;
end
for i=1:size(reg_1,1)
  if mod(i,1000)==0
    i
  end
  data_1(i).Sequence = reg_1(i,:);
  data_1(i).Header = 'ignore';
  data_1(i).Quality = char(abs('e')*ones(1,auxData.readLength));
end

system(['mkdir /home/csfaculty/shental/CS/BAC/EMIRGE/software/',fileName])
delete(['~/CS/BAC/EMIRGE/csmiller-EMIRGE-85bf31a/try/try_',fileName,'.fastq'])
fastqwrite(['~/CS/BAC/EMIRGE/csmiller-EMIRGE-85bf31a/try/try_',fileName,'.fastq'],data_1)


% write bact
clear c
c = struct;
for i=1:length(bact)
  c(i) = struct;
end

for i=1:length(bact)
  c(i).Sequence = Sequence_uni{bact(i)};
  c(i).Header = Header_uni{bact(i)};
  c(i).Quality = char(abs('e')*ones(1,length(c(i).Sequence)));
end
delete(['~/CS/BAC/EMIRGE/csmiller-EMIRGE-85bf31a/try/listOfBact_',fileName,'.fastq'])
fastqwrite(['~/CS/BAC/EMIRGE/csmiller-EMIRGE-85bf31a/try/listOfBact_',fileName,'.fastq'],c);
%%%%%%%%%%%%%%%%%%5




['~/CS/BAC/EMIRGE/csmiller-EMIRGE-85bf31a/emirge_amplicon.py /home/csfaculty/shental/CS/BAC/EMIRGE/software/',fileName,' -1 ~/CS/BAC/EMIRGE/csmiller-EMIRGE-85bf31a/try/try_',fileName,'.fastq  -f ~/CS/BAC/EMIRGE/csmiller-EMIRGE-85bf31a/SSU_candidate_db.fasta -b   ~/CS/BAC/EMIRGE/csmiller-EMIRGE-85bf31a/BWT/SSU_candidate_db_btindex -l 100 -i 1 -s 1']

['~/CS/BAC/EMIRGE/csmiller-EMIRGE-85bf31a/emirge_rename_fasta.py /home/csfaculty/shental/CS/BAC/EMIRGE/software/',fileName,'/iter.40 > /home/csfaculty/shental/CS/BAC/EMIRGE/software/',fileName,'/renamed.fasta  ']

  
