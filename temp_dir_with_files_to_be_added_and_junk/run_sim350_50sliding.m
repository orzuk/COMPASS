if 1==2 % creating reads

% sliding50 create reads
clear
createUniqueReads_sim350_50sliding('sim50_50silding_0p1perNoise')
createUniqueReads_sim350_50sliding('sim50_50silding_1perNoise')
createUniqueReads_sim350_50sliding('sim50_50silding_5perNoise')
createUniqueReads_sim350_50sliding('sim50_50silding_noNoise')

%%%%%%%%%%%%%%%%%%55
  
end

if 1==2 % prepare the files
  testSimFromReadsLikeRealData_50sliding_454({'sim50_50silding_0p1perNoise','sim50_50silding_1perNoise','sim50_50silding_5perNoise','sim50_50silding_noNoise'},'sliding50');
end

outputName = {'sim50_50silding_0p1perNoise','sim50_50silding_1perNoise','sim50_50silding_5perNoise','sim50_50silding_noNoise'};

for i=1:length(outputName)
  unix(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_',outputName{i},'/run_sol_',outputName{i},'_noCorrection_50'])
  pause(5)
end

for i=1:length(outputName)
  unix(['cp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/sol_',outputName{i},'/sol_',outputName{i},'_noCorrection_50/sol_',outputName{i},'_noCorrection_50.mat',' /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/sliding50_454/res_slising50_454'])

end


%%%%%%%%%%%%%%%%%%%%%5555
% compare
cd ~/CS/BAC/12Samples/454/idx/blast-2.2.17/bin

clear
userDir = getuserdir;
basicDir = [userDir,'/CS/BAC/sliding50_454/tmpCheck'];

outputName = {'sim350_50silding_0p1perNoise','sim350_50silding_1perNoise','sim350_50silding_5perNoise','sim350_50silding_noNoise'};

matlabpool local 6
for i=3%1:length(outputName) 
  seed = sum(100*clock);
  rand('seed',seed); 
  load(['~/CS/BAC/sliding50_454/',outputName{i}]); 
  clear A
  IND = 1:size(reads_uni,1);
  parfor k=1:size(reads_uni,1)
  
    currSeq = reads_uni(k,:);
    % changed for parallel
    % compare to without ambiguous
    A(k) = myFindSeqInDB_compareTo350WithoutAmbiguous(currSeq,basicDir,i,k);
  end
  bestMatch = A;clear A
  save(['~/CS/BAC/sliding50_454/data/reads_samples_',outputName{i},'_andMatch_basedDataBaseOf350WithoutAmbiguous'],'reads_uni','reads_uni_freq','bestMatch')
end
matlabpool close  


%%%%%%%%%%%%%%%%%%%%%%%%
% test cluster version

cd

clear
userDir = getuserdir;

  
  
basicDir = [userDir,'/CS/BAC/testAlign'];

basicName = 'testAlign1'; % should be a differnet name for each part
databaseName = [userDir,'/CS/BAC/12Samples/454/idx/old_fasta_files_without_ambiguous/454idx.fa'];

% use path
blastPath = [userDir,'/CS/BAC/12Samples/454/idx/blast-2.2.17/bin/'];
bestMatch = myFindSeqInDB_generic(reads_uni,basicName,basicDir,databaseName,blastPath,1,2);
  

  
auxData = struct;
auxData.userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
auxData.basicDir = [auxData.userDir,'/CS/BAC/testSmall'];
auxData.basicName = 'testAlign1';
auxData.databaseName = [auxData.userDir,'/CS/BAC/12Samples/454/idx/old_fasta_files_without_ambiguous/454idx.fa'];
auxData.blastPath = [auxData.userDir,'/CS/BAC/12Samples/454/idx/blast-2.2.17/bin/'];
auxData.dataDir = [auxData.userDir,'/CS/BAC/sliding50_454'];
auxData.dataName = ['smallData'];
auxData.numProcessors = 5;
auxData.queueName = 'hour';

unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/alignToDatabase/tmpRuns/alignTmp_',auxData.dataName])
unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/alignToDatabase/dataForSim/alignTmp_',auxData.dataName])
unix(['mkdir ',auxData.basicDir,';mkdir ',auxData.basicDir,'/tmp'])



alignToDatabse(auxData)


idx = myFindSeqInDB_generic(auxData.basicName,auxData.basicDir,auxData.databaseName,auxData.blastPath,1,2,saveFileName,auxData.dataDir,auxData.dataName)






%%%%%%%%%%%%%%
% read length 20

clear
readDir = 'sliding20_454';
readLength = 20;
createUniqueReads_sim350_sliding(['~/CS/BAC/',readDir],'sim20_20silding_noNoise',readLength)
createUniqueReads_sim350_sliding(['~/CS/BAC/',readDir],'sim20_20silding_0p1perNoise',readLength)
createUniqueReads_sim350_sliding(['~/CS/BAC/',readDir],'sim20_20silding_1perNoise',readLength)
createUniqueReads_sim350_sliding(['~/CS/BAC/',readDir],'sim20_20silding_5perNoise',readLength)


if 1==2 % prepare the files
  outputName = {'sim20_20silding_noNoise','sim20_20silding_0p1perNoise','sim20_20silding_1perNoise','sim20_20silding_5perNoise'};
  testSimFromReadsLikeRealData_sliding_454(outputName,readDir,20);
end


for i=1:length(outputName)
  unix(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_',outputName{i},'/run_sol_',outputName{i},'_noCorrection_',num2str(readLength)])
  pause(5)
end

unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/',readDir,'/res_slising',num2str(readLength),'_454'])
for i=1:length(outputName)
  unix(['cp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/sol_',outputName{i},'/sol_',outputName{i},'_noCorrection_',num2str(readLength),'/sol_',outputName{i},'_noCorrection_',num2str(readLength),'.mat',' /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/',readDir,'/res_slising',num2str(readLength),'_454'])

end

