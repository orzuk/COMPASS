%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% run with REVERSE

%testData12SamplesSolexa_WITH_REVERSE({'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S3','S4'},100);testData12SamplesSolexa_WITH_REVERSE({'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S3','S4'},50)

% no correction
currFileName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/run_no_correction_WITH_REVERSE.txt';
fid = fopen(currFileName,'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
fprintf(fid,'mkdir $TMPDIR\n');
fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
fprintf(fid,'\n');  

sampleName = {'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S3','S4'};
for i=1:length(sampleName)
  fprintf(fid,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_WITH_REVERSE_%s/run_sol_WITH_REVERSE_%s_noCorrection_100\n',sampleName{i},sampleName{i});
end
for i=1:length(sampleName)
  fprintf(fid,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_WITH_REVERSE_%s/run_sol_WITH_REVERSE_%s_noCorrection_50\n',sampleName{i},sampleName{i});
end

fprintf(fid,'rm -fr $TMPDIR\n');
fclose(fid);
unix(['chmod 0700 ',currFileName]);



% with correction
currFileName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/run_with_correction_WITH_REVERSE.txt';
fid = fopen(currFileName,'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
fprintf(fid,'mkdir $TMPDIR\n');
fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
fprintf(fid,'\n');  

sampleName = {'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S3','S4'};
for i=1:length(sampleName)
  fprintf(fid,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_WITH_REVERSE_%s/run_sol_WITH_REVERSE_%s_withCorrection_100\n',sampleName{i},sampleName{i});
end
for i=1:length(sampleName)
  fprintf(fid,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_WITH_REVERSE_%s/run_sol_WITH_REVERSE_%s_withCorrection_50\n',sampleName{i},sampleName{i});
end

fprintf(fid,'rm -fr $TMPDIR\n');
fclose(fid);
unix(['chmod 0700 ',currFileName]);

% find sol_WITH_REVERSE_* -name "*.mat" >> list_WITH_REVERSE
% basj
% for i in `cat list_WITH_REVERSE`; do cp $i resSolexa_WITH_REVERSE;done

1
pause

%%%%%%%%%%%%%%%%%%%%%555
% find Wolbachia

clear
load ~/CS/BAC/primers750/bac16s_primers750


%load /home/csfaculty/shental/CS/BAC/resNoCorrection/sol_WITH_REVERSE_M4_noCorrection_100
%load ~/CS/BAC/12Samples/Solexa/data/readsSolexa_WITH_REVERSE_100_sample_M4.mat

%load /home/csfaculty/shental/CS/BAC/resNoCorrection/sol_WITH_REVERSE_O7_noCorrection_100
%load ~/CS/BAC/12Samples/Solexa/data/readsSolexa_WITH_REVERSE_100_sample_O7.mat

load /home/csfaculty/shental/CS/BAC/resNoCorrection/sol_WITH_REVERSE_S4_noCorrection_100
load ~/CS/BAC/12Samples/Solexa/data/readsSolexa_WITH_REVERSE_100_sample_S4.mat


userDir = getuserdir;
basicSeqNameDir = [userDir,'/CS/BAC/primers750/datNoNonACGT/packed64/'];
basicSeqKey= [userDir,'/CS/BAC/primers750/datNoNonACGT/keyNoNonACGT_primers750'];

tmpInd = store_kp{end};
%tmpInd = find(abs(found{end})>10^-3);

[junk,ind] = sort(abs(found{end}),'descend');
a = find(junk>10^-3);
junk= junk(a);
ind = ind(a);

Header_750{ind(1)}


tmpInd = ind;

readLength = 100;
[normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(readLength,tmpInd,basicSeqNameDir,basicSeqKey);


dataIn = struct;
[fracRelevantReads,sumRelevantReads(i)] = currReadsFourth(uniqueReads,uniqueReads_length,values,0,dataIn);



fracRelevantReads = myCorrectReadsNew3_for_SL4(values,uniqueReads,uniqueReads_length');

a = sum(normalizedBac,1); 

b = find(fracRelevantReads);
sum(fracRelevantReads(b))

% add - number of reads - add in code - save this
% 
%%%%%%%5

%%%%%%%%%%%%%%%%%55
% compare olutions

userDir = getuserdir;
dirName = '~/CS/BAC/resSolexa_WITH_REVERSE/'
auxData.basicSeqNameDir = [userDir,'/CS/BAC/primers750/datNoNonACGT/packed64/'];
auxData.basicSeqKey = [userDir,'/CS/BAC/primers750/datNoNonACGT/keyNoNonACGT_primers750'];
auxData.currDir = '~/CS/BAC/junk/';
auxData.userDir = userDir;
CreateMothurDistTwoSolutions([dirName,'sol_WITH_REVERSE_S7_withCorrection_50.mat'],[dirName,'sol_WITH_REVERSE_S7_noCorrection_50.mat'],auxData);


CreateMothurDistTwoSolutions([dirName,'sol_WITH_REVERSE_S1_noCorrection_50.mat'],[dirName,'sol_WITH_REVERSE_S1_withCorrection_100.mat'],auxData);


compareTwoResults([dirName,'sol_WITH_REVERSE_S1_noCorrection_50.mat'],[dirName,'sol_WITH_REVERSE_S1_withCorrection_100.mat'],0.01,auxData)

compareTwoResults([dirName,'sol_WITH_REVERSE_O10_noCorrection_100.mat'],[dirName,'sol_WITH_REVERSE_O10_withCorrection_100.mat'],0.01,auxData)

