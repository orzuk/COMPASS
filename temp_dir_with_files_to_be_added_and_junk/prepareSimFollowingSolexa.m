%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
load ~/CS/BAC/simFollowingSolexa/reads-mix1

% find uniqueReads

r = int2nt(unpack_seqs(reads,100,64));
readChar = cell(size(r,1),1);
for i=1:size(r,1)
  readChar{i} = r(i,:);
end

m4sort=sort(readChar);
[uni_reads,ia]=unique(m4sort,'first');
[~,ib]=unique(m4sort,'last');
uniqueReads_length = ib-ia+1;

uniqueReads = pack_seqs(uni_reads,64);
uniqueReads = cell2mat(uniqueReads);

save ~/CS/BAC/simFollowingSolexa/data_reads_mix1  uniqueReads_length uniqueReads
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare the reads - all - run with correction


%two runs
datasetName = 'primers750_primer_tail';
testSimFromReadsLikeRealData('test_mix_100_1','simFollowingSolexa/data_reads_mix1',datasetName)

%%%%
clear
load ~/CS/BAC/sol_test_mix_100_1_test_mix_100_1_noCorrection_100
load ~/CS/BAC/simFollowingSolexa/reads-mix1
load ~/CS/BAC/simFollowingSolexa/data_reads_mix1
readLength = 100;

find(abs(found{3})>10^-3)

userDir = getuserdir;
datasetName = 'primers750_primer_tail';
basicSeqNameDir = [userDir,'/CS/BAC/',datasetName,'/datNoNonACGT/packed64/'];
basicSeqKey = [userDir,'/CS/BAC/',datasetName,'/datNoNonACGT/keyNoNonACGT_',datasetName];

tmpInd = find(abs(found{3})>10^-3);
[normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(readLength,tmpInd,basicSeqNameDir,basicSeqKey);
  
  
dataIn = struct;
[fracRelevantReads,sumRelevantReads] = currReadsFourth(uniqueReads,uniqueReads_length,values,0,dataIn);

[x1_2_nor] = testL1_2(normalizedBac./max(fracRelevantReads),fracRelevantReads./max(fracRelevantReads));    
[x1_2_no_nor] = testL1_2(normalizedBac,fracRelevantReads);    



%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
load ~/CS/BAC/simFollowingSolexa/reads-mix100

% find uniqueReads

r = int2nt(unpack_seqs(reads,100,64));
readChar = cell(size(r,1),1);
for i=1:size(r,1)
  readChar{i} = r(i,:);
end

m4sort=sort(readChar);
[uni_reads,ia]=unique(m4sort,'first');
[~,ib]=unique(m4sort,'last');
uniqueReads_length = ib-ia+1;

uniqueReads = pack_seqs(uni_reads,64);
uniqueReads = cell2mat(uniqueReads);

save ~/CS/BAC/simFollowingSolexa/data_reads_mix100  uniqueReads_length uniqueReads
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
load ~/CS/BAC/sol_test_mix_100_test_mix_100_noCorrection_100
load ~/CS/BAC/simFollowingSolexa/reads-mix100

userDir = getuserdir;
datasetName = 'primers750_primer_tail';
basicSeqNameDir = [userDir,'/CS/BAC/',datasetName,'/datNoNonACGT/packed64/'];
basicSeqKey = [userDir,'/CS/BAC/',datasetName,'/datNoNonACGT/keyNoNonACGT_',datasetName];

tmpInd = find(abs(found{end})>10^-3);
[normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(readLength,tmpInd,basicSeqNameDir,basicSeqKey);

 
dataIn = struct;
[fracRelevantReads,sumRelevantReads] = currReadsFourth(uniqueReads,uniqueReads_length,values,0,dataIn);

[x1_2_nor] = testL1_2(normalizedBac./max(fracRelevantReads),fracRelevantReads./max(fracRelevantReads));    
[x1_2_no_nor] = testL1_2(normalizedBac,fracRelevantReads);    

f = zeros(212040,1);
f(tmpInd) = x1_2_nor;

f2 = zeros(212040,1);
f2(tmpInd) = abs(found{end}(tmpInd));

c = zeros(212040,1);
c(bact) = freq;



% prepare the reads - all - run with correction

dist
matrix


%two runs
datasetName = 'primers750_primer_tail';
testSimFromReadsLikeRealData('test_mix_100','simFollowingSolexa/data_reads_mix100',datasetName)

%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%

% 3 rand
clear

fileName = 'data_reads_3rand';
inputName = 'reads-mix3-rand';
outputName = 'test_mix_3rand';


createUniqueReadsSimFollowingSolexa(['~/CS/BAC/simFollowingSolexa/',inputName],['~/CS/BAC/simFollowingSolexa/',fileName])

% copy to br

datasetName = 'primers750_primer_tail';
testSimFromReadsLikeRealData(outputName,['simFollowingSolexa/',fileName],datasetName)

load(['~/CS/BAC/sol_',outputName,'_',outputName,'_noCorrection_100']); % 
load(['~/CS/BAC/simFollowingSolexa/',fileName]) 
load(['~/CS/BAC/simFollowingSolexa/',inputName])
find(abs(found{end})>10^-3)

userDir = getuserdir;
datasetName = 'primers750_primer_tail';
basicSeqNameDir = [userDir,'/CS/BAC/',datasetName,'/datNoNonACGT/packed64/'];
basicSeqKey = [userDir,'/CS/BAC/',datasetName,'/datNoNonACGT/keyNoNonACGT_',datasetName];
readLength = 100;

tmpInd = find(abs(found{end})>10^-3);
[normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(readLength,tmpInd,basicSeqNameDir,basicSeqKey);

 
dataIn = struct;
[fracRelevantReads,sumRelevantReads] = currReadsFourth(uniqueReads,uniqueReads_length,values,0,dataIn);

[x1_2_nor] = testL1_2(normalizedBac./max(fracRelevantReads),fracRelevantReads./max(fracRelevantReads));    
[x1_2_no_nor] = testL1_2(normalizedBac,fracRelevantReads);    

f = zeros(212040,1);
f(tmpInd) = x1_2_nor;

f2 = zeros(212040,1);
f2(tmpInd) = abs(found{end}(tmpInd));

c = zeros(212040,1);
c(bact) = freq;

clf
plot(c,f,'b.')
hold on
plot(c,f2,'kx')
legend('l1','l2')
title(inputName)
print('-dpdf',['~/CS/BAC/',inputName])  


figure(2)
plot(c,f,'kx')
line([0,max(c)],[0 max(c)])
xlabel('correct');ylabel('l1')
figure(3)
plot(c,f2,'kx')
line([0,max(c)],[0 max(c)])
xlabel('correct');ylabel('l2')
figure(4)
plot(f,f2,'.')
line([0,max(c)],[0 max(c)])
xlabel('l1');ylabel('l2')

save(['~/CS/BAC/simFollowingSolexa/res_',inputName],'f','f2','c')


%%%%%%
%%%%%
% 100-uni

clear
fileName = 'data_reads_100uni';
inputName = 'reads-mix100-uni';
outputName = 'test_mix_100uni';

createUniqueReadsSimFollowingSolexa(['~/CS/BAC/simFollowingSolexa/',inputName],['~/CS/BAC/simFollowingSolexa/',fileName])

% copy to br

datasetName = 'primers750_primer_tail';
testSimFromReadsLikeRealData(outputName,['simFollowingSolexa/',fileName],datasetName)

load(['~/CS/BAC/sol_',outputName,'_',outputName,'_noCorrection_100']); % 
load(['~/CS/BAC/simFollowingSolexa/',fileName]) 
load(['~/CS/BAC/simFollowingSolexa/',inputName])
find(abs(found{end})>10^-3)


userDir = getuserdir;
datasetName = 'primers750_primer_tail';
basicSeqNameDir = [userDir,'/CS/BAC/',datasetName,'/datNoNonACGT/packed64/'];
basicSeqKey = [userDir,'/CS/BAC/',datasetName,'/datNoNonACGT/keyNoNonACGT_',datasetName];
readLength = 100;

tmpInd = find(abs(found{end})>10^-3);
[normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(readLength,tmpInd,basicSeqNameDir,basicSeqKey);

 
dataIn = struct;
[fracRelevantReads,sumRelevantReads] = currReadsFourth(uniqueReads,uniqueReads_length,values,0,dataIn);

[x1_2_nor] = testL1_2(normalizedBac./max(fracRelevantReads),fracRelevantReads./max(fracRelevantReads));    
[x1_2_no_nor] = testL1_2(normalizedBac,fracRelevantReads);    

f = zeros(212040,1);
f(tmpInd) = x1_2_nor;

f2 = zeros(212040,1);
f2(tmpInd) = abs(found{end}(tmpInd));

c = zeros(212040,1);
c(bact) = freq;

clf
plot(c,f,'b.')
hold on
plot(c,f2,'kx')
legend('l1','l2')
title(inputName)
print('-dpdf',['~/CS/BAC/',inputName])  


figure(2)
plot(c,f,'kx')
line([0,max(c)],[0 max(c)])
xlabel('correct');ylabel('l1')
figure(3)
plot(c,f2,'kx')
line([0,max(c)],[0 max(c)])
xlabel('correct');ylabel('l2')
figure(4)
plot(f,f2,'.')
line([0,max(c)],[0 max(c)])
xlabel('l1');ylabel('l2')


save(['~/CS/BAC/simFollowingSolexa/res_',inputName],'f','f2','c')













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%
%%%%%%
%%%%%
% 24April


clear
basicName = {'reads-S1-noise-0','reads-S1-noise-100','reads-S1-noise-10','reads-S1-noise-1','reads-S1-noise-S1'};
auxDataTest = struct;

%auxDataTest.datasetName = 'primers750_primer_tail';
%auxDataTest.numBACtoConsider = 212040;
%auxDataTest.dataSetPrimers750Flag = 1;

auxDataTest.datasetName = 'fullDatabse';
auxDataTest.numBACtoConsider = 410849;
auxDataTest.dataSetPrimers750Flag = 0;


auxDataTest.readLength = 100;
auxDataTest.currNumProcessors = 40;
auxDataTest.basicDirNameForFile = 'sim24April';

clear fileName inputName outputName
for i=1:length(basicName)
  fileName{i} = ['myFormat_',basicName{i}];
  inputName{i} = basicName{i};
  outputName{i} = ['test_',basicName{i}];
end

for i=1:length(basicName)
  createUniqueReadsSimFollowingSolexa(['~/CS/BAC/simFollowingSolexa/',auxDataTest.basicDirNameForFile,'/',inputName{i}],['~/CS/BAC/simFollowingSolexa/',auxDataTest.basicDirNameForFile,'/',fileName{i}],auxDataTest.readLength)
end


% copy to br
PWD = pwd;
cd(['~/CS/BAC/simFollowingSolexa/',auxDataTest.basicDirNameForFile])
clear w
w = 'scp '
for i=1:length(basicName)
  w = [w,' ',fileName{i},'.mat '];
end
w = [w,' orzuk@tin.broadinstitute.org:/seq/orzuk2/compressed_sensing/metagenomics/next_gen','/CS/BAC/simFollowingSolexa/',auxDataTest.basicDirNameForFile,'/'];

cd(PWD);


% prepare the files

for i=1:length(outputName)
  testSimFromReadsLikeRealData_generic(outputName{i},['simFollowingSolexa/sim24April/',fileName{i}],auxDataTest)
end

%%%%%5
% run the files

for i=1:length(outputName)
  unix(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_',outputName{i},'_',outputName{i},'/run_sol_',outputName{i},'_',outputName{i},'_noCorrection_',num2str(auxDataTest.readLength)])
end

% collect results
for i=1:length(outputName)
  unix(['cp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/sol_',outputName{i},'_',outputName{i},'/sol_',outputName{i},'_',outputName{i},'_noCorrection_',num2str(auxDataTest.readLength),'/sol_',outputName{i},'_',outputName{i},'_noCorrection_',num2str(auxDataTest.readLength),'.mat',' /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/simFollowingSolexa/',auxDataTest.basicDirNameForFile,'/results/'])
end

% copy to sol
unix(['cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/simFollowingSolexa/',auxDataTest.basicDirNameForFile,';tar cf results_', ...
      auxDataTest.basicDirNameForFile,'.tar ','/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/simFollowingSolexa/',auxDataTest.basicDirNameForFile,'/results/'])

% copy locally
unix(['scp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/simFollowingSolexa/',auxDataTest.basicDirNameForFile,'/results_', ...
      auxDataTest.basicDirNameForFile,'.tar ',...
      'shental@sol.cslab.openu.ac.il:~/CS/BAC/simFollowingSolexa/',auxDataTest.basicDirNameForFile])
% run locally

for i=1:length(outputName)
  extractResSim_generic(outputName{i},inputName{i},fileName{i},auxDataTest)
  drawnow
  pause
end



