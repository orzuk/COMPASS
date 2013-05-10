% 10.5 solves the forward and reverse mixup - add cases in which the reads of forward and reverse are similar

% create reads
clear

dirName = 'opt3_MergeForwardReverse';
typeName = '_unifreq_corr_opt3_th50';
datasetName = 'primers750_primer_tail';
basicName = {'S1','S2','S3','S4','M1','M2','M3','M4','S7','S10','O7','O10'};

readLength = 100;
thresholdForReverse = 0; % POS>0 for forward and POS<0 for reverse
for i=1:length(basicName)
  createUniqueReads_12SamplesSolexa_eq_unicov_MergeForwardReverse(dirName,readLength,basicName{i},typeName,thresholdForReverse)
end


%%%%%%%%%%%%
% create the files for Amnon - without the opt3 data
for i=1:length(basicName)
  clear uniqueReads uni_read_seq freq_ForwardAndReverse
  load(['~/CS/BAC/12Samples/Solexa/data/',dirName,'/data_',dirName,'_',basicName{i}],'uniqueReads','uni_read_seq','freq_ForwardAndReverse');
  save(['~/CS/BAC/12Samples/Solexa/data/',dirName,'/data_NO_OPT3_',dirName,'_',basicName{i}],'uniqueReads','uni_read_seq','freq_ForwardAndReverse');
end


% prepare the run files equalization

auxDataTest = struct;
auxDataTest.datasetName = 'primers750_primer_tail';
auxDataTest.numBACtoConsider = 212040;
auxDataTest.dataSetPrimers750Flag = 1;

auxDataTest.readLength = 100;

auxDataTest.currNumProcessors = 20;
auxDataTest.basicDirNameForFile = 'opt3_MergeForwardReverse';

clear fileName inputName outputName
for i=1:length(basicName)
  fileName{i} = ['data_',auxDataTest.basicDirNameForFile,'_',basicName{i}];
  outputName{i} = ['test_',auxDataTest.basicDirNameForFile,'_',basicName{i}];
end

% prepare the files
for i=1:length(outputName)
  testSimFromReadsLikeRealData_generic2(outputName{i},['12Samples/Solexa/data/',auxDataTest.basicDirNameForFile,'/',fileName{i}],auxDataTest)
end

%%%%%5
% run the files
for i=1:length(outputName)
  unix(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_',outputName{i},'/run_sol_',outputName{i},'_noCorrection_',num2str(auxDataTest.readLength)])
end

% collect results
unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/12Samples/Solexa/results/',auxDataTest.basicDirNameForFile])
for i=1:length(outputName)
  unix(['cp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/sol_',outputName{i},'/sol_',outputName{i},'_noCorrection_',num2str(auxDataTest.readLength),'/sol_',outputName{i},'_noCorrection_',num2str(auxDataTest.readLength),'.mat',' /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/12Samples/Solexa/results/',auxDataTest.basicDirNameForFile])
end

unix(['cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/12Samples/Solexa/results/',';tar cf results_', ...
      auxDataTest.basicDirNameForFile,'.tar ',auxDataTest.basicDirNameForFile])

['scp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/12Samples/Solexa/results/results_',auxDataTest.basicDirNameForFile,'.tar shental@sol.cslab.openu.ac.il:~/CS/BAC/12Samples/Solexa/results/']
% run locally

resultsDir = ['~/CS/BAC/12Samples/Solexa/results'];
for i=1:length(outputName)
  readDataFileName = [userDir,'/CS/BAC/12Sampreles/Solexa/data/',auxDataTest.basicDirNameForFile,'/',fileName{i}]
  extractResRealData_generic(outputName{i},readDataFileName,resultsDir,auxDataTest)
  drawnow
  %pause
end

for i=1:length(outputName)
  plotResRealData_generic(resultsDir,outputName{i},auxDataTest)
  drawnow
  pause
end


i=4
res = load([resultsDir,'/',auxDataTest.basicDirNameForFile,'/res_',outputName{i}]);

load ~/CS/BAC/primers750_primer_tail/bac16s_primers750_primer_tail_full_without_ambiguous

a = find(res.x1_2_nor>10^-3);
[junk,i1] = sort(res.x1_2_nor(a),'descend');


for i=1:length(i1)
  for j=1:length(Header_750_tail{res.tmpInd(a(i1(i)))})
    Header_750_tail{res.tmpInd(a(i1(i)))}{j}
    junk(i)
  end
  disp('%%%%%')
  pause
end


basicSeqNameDir = [userDir,'/CS/BAC/',auxDataTest.datasetName,'/datNoNonACGT/packed64/'];
basicSeqKey= [userDir,'/CS/BAC/',auxDataTest.datasetName,'/datNoNonACGT/keyNoNonACGT_',auxDataTest.datasetName];
load(basicSeqKey)

[HeaderAll,SequenceAll] = loadSeq(res.tmpInd,basicSeqNameDir,basicSeqKey);
ham = zeros*ones(length(res.tmpInd));

for i=1:length(res.tmpInd)-1
  i
  seq_tmpInd = int2nt(unpack_seqs(SequenceAll{i},len_uni(res.tmpInd(i)),64));
  for j=i+1:length(res.tmpInd)
    seq_rec = int2nt(unpack_seqs(SequenceAll{j},len_uni(res.tmpInd(j)),64));
    [cscore,algn]=swalign(seq_tmpInd,seq_rec,'Alphabet','NT');
    ham(i,j) = length(find(algn=='|'));
    
    %ham(i,j)
    %pause
  end
  length(find((ham)))
end

ham = ham+ham';

for i=1:length(res.tmpInd)
  ham(i,i) = len_uni(res.tmpInd(i));
end

currH = cell(1,length(res.tmpInd));
for i=1:length(res.tmpInd)
  b = Header_750_tail{res.tmpInd(i)}{1}
  currH{i} = b(20:min(60,length(b)));
end

figure(1);clf
imagesc(ham)
set(gca,'xtick',1:length(currH),'xticklabel',currH,'position',[0.1300    0.2100    0.7750    0.7150 ],'ytick',1:length(currH))
rotateticklabel(gca,45)

for i=1:length(i1)
  
  figure(2);clf
  plot(ham(a(i1(i)),a(i1(i)))-ham(a(i1(i)),:),'.');
  hold on
  plot(a(i1(i)),0,'ro')
  Header_750_tail{res.tmpInd(a(i1(i)))}{1}
  set(gca,'xtick',1:length(currH),'xticklabel',currH,'position',[0.1300    0.2100    0.7750    0.7150 ],'fontweight','bold','fontsize',20)
  rotateticklabel(gca,30);
  set(gca,'fontweight','bold','fontsize',20)
  pause
end




a = find(res.f_l2>10^-3);
[junk,i1] = sort(res.f_l2(a),'descend');
for i=1:length(i1)
  for j=1:length(Header_750_tail{a(i1(i))})
    Header_750_tail{a(i1(i))}{j}
  end
  disp('%%%%%')
  pause
end

