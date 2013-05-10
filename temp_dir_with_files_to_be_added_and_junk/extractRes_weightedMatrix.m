function extractRes_weightedMatrix(outputName,readDataFileName,profileName,resultsDir,auxDataTest,filterFlag)
%keyboard

load([resultsDir,'/',auxDataTest.basicDirNameForFile,'/sol_',outputName,'_noCorrection_',num2str(auxDataTest.readLength)],'found'); 
load([readDataFileName],'uniqueReads','uniqueReads_length_forward','uniqueReads_length_reverse')
load(profileName)

uniqueReads_length = uniqueReads_length_forward+uniqueReads_length_reverse;


find(abs(found{end})>10^-3)

userDir = getuserdir;

if strcmp(auxDataTest.datasetName,'fullDatabse')
  basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
  basicSeqKey= [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
else
  basicSeqNameDir = [userDir,'/CS/BAC/',auxDataTest.datasetName,'/datNoNonACGT/packed64/'];
  basicSeqKey= [userDir,'/CS/BAC/',auxDataTest.datasetName,'/datNoNonACGT/keyNoNonACGT_',auxDataTest.datasetName];
  
end

tmpInd = find(abs(found{end})>10^-3);

%keyboard
[normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(auxDataTest.readLength,tmpInd,basicSeqNameDir,basicSeqKey);

dataIn = struct;
[fracRelevantReads,sumRelevantReads] = currReadsFourth(uniqueReads,uniqueReads_length,values,0,dataIn);
%keyboard 

numBACtoConsider = auxDataTest.numBACtoConsider;
%keyboard

if filterFlag==1
  auxDataTest.fvalFor = fvalFor;
  auxDataTest.fvalRev = fvalRev;
  disp('using filtered noise')
else
  auxDataTest.fvalFor = valFor;
  auxDataTest.fvalRev = valRev;
  disp('using non-filtered noise')
end


matrix = weightedMatrix(auxDataTest,tmpInd,basicSeqNameDir,basicSeqKey);

ratio = max(sum(matrix,1)./max(fracRelevantReads));
tic
freqWeightedMatrix = testL1_2(matrix./max(fracRelevantReads)/ratio,fracRelevantReads./max(fracRelevantReads)/ratio);
toc
%keyboard

%tic
%[x{i}]=runOneGroupOf1000ForCompilationFourth(normalizedBac,fracRelevantReads);    
%toc



save([resultsDir,'/',auxDataTest.basicDirNameForFile,'/resWeightedMatrixBased_',outputName],'fracRelevantReads','normalizedBac','matrix','values','freqWeightedMatrix','tmpInd')

plot(freqWeightedMatrix,'.')
