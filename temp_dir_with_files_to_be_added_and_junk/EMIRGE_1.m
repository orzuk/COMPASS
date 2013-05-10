clear

numBac = 410849;
bact = ceil(rand(1,20)*numBac);
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
for i=bact
  
  a = find(Header_uni{i}==' ');
  Header_uni{i}(a(1)+1:a(2)-1)
  %pause
end

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

delete('~/CS/BAC/EMIRGE/csmiller-EMIRGE-85bf31a/try/try_1.fastq')
fastqwrite('~/CS/BAC/EMIRGE/csmiller-EMIRGE-85bf31a/try/try_1.fastq',data_1)
%fastqwrite('~/CS/BAC/EMIRGE/csmiller-EMIRGE-85bf31a/try/try_2.fastq',data_2)


./emirge_amplicon.py /home/csfaculty/shental/CS/BAC/EMIRGE/software/t1 -1 ~/CS/BAC/EMIRGE/csmiller-EMIRGE-85bf31a/try/try_1.fastq  -f ~/CS/BAC/EMIRGE/csmiller-EMIRGE-85bf31a/SSU_candidate_db.fasta -b   ~/CS/BAC/EMIRGE/csmiller-EMIRGE-85bf31a/BWT/SSU_candidate_db_btindex -l 100 -i 1 -s 1

./emirge_rename_fasta.py /home/csfaculty/shental/CS/BAC/EMIRGE/software/t1/iter.40 > /home/csfaculty/shental/CS/BAC/EMIRGE/software/t1/renamed.fasta  
  
