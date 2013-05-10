clear
load ~/CS/BAC/primers750_primer_tail/bac16s_primers750_primer_tail_full_without_ambiguous

fileName = '~/CS/BAC/primers750_primer_tail/fastaFile/bac16s_primers750_primer_tail_full_without_ambiguous.fa'
delete([fileName])

for i=1:length(Header_750_tail)
  fastawrite(fileName,num2str(i),Sequence_750_tail{i});
end


blastformat('Inputdb','~/CS/BAC/primers750_primer_tail/fastaFile/bac16s_primers750_primer_tail_full_without_ambiguous.fa','Protein','false');


