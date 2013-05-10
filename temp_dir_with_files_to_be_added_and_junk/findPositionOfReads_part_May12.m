function findPositionOfReads_part_May12(sampleName,readLength,saveFileName,indexFirst,indexLast)
% rmpath(genpath('/home/csfaculty/shental/Weizmann/ph2users/fenoam/'));rmpath('/home/csfaculty/shental/matlab');cd ~/CS/mFilesBAC;addpath('~/CS/mFilesBAC');mcc -v -m findPositionOfReads_part_May12.m -d /homes/csfaculty/shental/CS/BAC/cvx -a get_mean_posalign_v3

disp('run get_mean_posalign_v3')

if ~isdeployed
  userDir = getuserdir
else
  userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
  readLength = str2num(readLength);
  indexFirst = str2num(indexFirst);
  indexLast = str2num(indexLast);

end


load([userDir,'/CS/BAC/12Samples/Solexa/data/forAlign/primer750_primer_tail_charVector_afterIncludingRightTail']) 

data = load([userDir,'/CS/BAC/12Samples/Solexa/data/opt3_th50_mar25_2012/illumina_reads',num2str(readLength),'_',sampleName,'_unifreq_corr_opt3_th50'],'uni_read_seq','POS')


for i=indexFirst:indexLast
  tmp_read = data.uni_read_seq(i,:); % RELEVANT ONLY FOR THE CURRENT ONE IN WHICH ALL ARE FORWARD
  if data.POS(i)<0
    tmp_read = seqrcomplement(tmp_read);
  end
  [pos(i) rev_pos(i) sd(i) rev_sd(i)]=get_mean_posalign_v3(tmp_read,database1,database_size,bacseq_len);
end

save([saveFileName],'pos','rev_pos', 'sd', 'rev_sd','indexFirst','indexLast')
