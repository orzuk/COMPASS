% testData12SamplesSolexa_br({'O7','O10','S7','S10'},100);testData12SamplesSolexa_br({'O7','O10','S7','S10'},50);testData12SamplesSolexa_br({'M1','M2','M3','M4'},100);testData12SamplesSolexa_br({'M1','M2','M3','M4'},50);testData12SamplesSolexa_br({'S1','S2','S4'},100);testData12SamplesSolexa_br({'S1','S2','S4'},50);

currFileName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/run_no_correction.txt';
fid = fopen(currFileName,'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
fprintf(fid,'mkdir $TMPDIR\n');
fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
fprintf(fid,'\n');  

sampleName = {'O7','O10','S7','S10'};
for i=1:length(sampleName)
  fprintf(fid,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_%s/run_sol_%s_noCorrection_100\n',sampleName{i},sampleName{i});
  fprintf(fid,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_%s/run_sol_%s_noCorrection_50\n',sampleName{i},sampleName{i});
end
fprintf(fid,'rm -fr $TMPDIR\n');
fclose(fid);
unix(['chmod 0700 ',currFileName]);

%%%%%%%%%%%%%%
% run the rest for no correction
currFileName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/run_no_correction2.txt';
fid = fopen(currFileName,'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
fprintf(fid,'mkdir $TMPDIR\n');
fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
fprintf(fid,'\n');  

sampleName = {'M1','M2','M3','M4','S1','S2','S4'};
for i=1:length(sampleName)
  fprintf(fid,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_%s/run_sol_%s_noCorrection_100\n',sampleName{i},sampleName{i});
  fprintf(fid,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_%s/run_sol_%s_noCorrection_50\n',sampleName{i},sampleName{i});
end
fprintf(fid,'rm -fr $TMPDIR\n');
fclose(fid);
unix(['chmod 0700 ',currFileName]);



% /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_O7/run_sol_O7_noCorrection_100



% find sol_* -name "*.mat" >> list
% basj
% for i in `cat list`; do cp $i temp;done


%%%%%%%%%%%%%%%%%%%555
% with correction

currFileName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/run_with_correction.txt';
fid = fopen(currFileName,'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
fprintf(fid,'mkdir $TMPDIR\n');
fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
fprintf(fid,'\n');  

sampleName = {'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S4'};
for i=1:length(sampleName)
  fprintf(fid,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_%s/run_sol_%s_withCorection_100\n',sampleName{i},sampleName{i});
end
for i=1:length(sampleName)
  fprintf(fid,'/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_%s/run_sol_%s_withCorection_50\n',sampleName{i},sampleName{i});
end

fprintf(fid,'rm -fr $TMPDIR\n');
fclose(fid);
unix(['chmod 0700 ',currFileName]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% run with REVERSE

%testData12SamplesSolexa_WITH_REVERSE({'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S4'},100);testData12SamplesSolexa_WITH_REVERSE({'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S4'},50)

% no correction
currFileName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/run_no_correction_WITH_REVERSE.txt';
fid = fopen(currFileName,'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,'TMPDIR=/tmp/${RANDOM}_{RANDOM}\n');
fprintf(fid,'mkdir $TMPDIR\n');
fprintf(fid,'export MCR_CACHE_ROOT=$TMPDIR\n');
fprintf(fid,'\n');  

sampleName = {'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S4'};
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

sampleName = {'O7','O10','S7','S10','M1','M2','M3','M4','S1','S2','S4'};
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

%%%%%%%%%%%%%%%%%%%%%555
% find Wolbachia

clear
load /home/csfaculty/shental/CS/BAC/resSolexa_WITH_REVERSE/sol_WITH_REVERSE_M4_withCorrection_100

load ~/CS/BAC/primers750/bac16s_primers750
load ~/CS/BAC/12Samples/Solexa/data/readsSolexa_WITH_REVERSE_100_sample_M4.mat


userDir = getuserdir;
basicSeqNameDir = [userDir,'/CS/BAC/primers750/datNoNonACGT/packed64/'];
basicSeqKey= [userDir,'/CS/BAC/primers750/datNoNonACGT/keyNoNonACGT_primers750'];

tmpInd = store_kp{2};
%tmpInd = find(abs(found{end})>10^-3);

[junk,ind] = sort(abs(found{end}),'descend');
a = find(junk>10^-3);
junk= junk(a);
ind = ind(a);

tmpInd = ind;

readLength = 100;
[normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(readLength,tmpInd,basicSeqNameDir,basicSeqKey);


dataIn = struct;
[fracRelevantReads,sumRelevantReads(i)] = currReadsFourth(uniqueReads,uniqueReads_length,values,0,dataIn);



fracRelevantReads = myCorrectReadsNew3_for_SL4(values,uniqueReads,uniqueReads_length');

a = sum(normalizedBac,1); 

b = find(fracRelevantReads);
sum(fracRelevantReads(b))
