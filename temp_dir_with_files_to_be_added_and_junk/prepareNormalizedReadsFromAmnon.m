% create files from Amnon's normalized read

% Amnon's file Len90.zip: contaims files NormRead* which are renamed Window90Norm_
clear
basicName = {'S1','S2','S3','S4','M1','M2','M3','M4','S7','S10','O7','O10'};
for i=2:length(basicName)
  unix(['mv /home/csfaculty/shental/CS/BAC/12Samples/Solexa/data/window_90_withNor_MergeForwardAndReverse/NormReads',basicName{i},'.mat /home/csfaculty/shental/CS/BAC/12Samples/Solexa/data/window_90_withNor_MergeForwardAndReverse/Window90Norm_',basicName{i},'.mat']);
end


%%%%% for normReadsf
clear
basicName = {'S1','S2','S3','S4','M1','M2','M3','M4','S7','S10','O7','O10'};

auxDataTest.datasetName = 'primers750_primer_tail';
auxDataTest.numBACtoConsider = 212040;
auxDataTest.dataSetPrimers750Flag = 1;

auxDataTest.readLength = 90;

auxDataTest.currNumProcessors = 20;
auxDataTest.basicDirNameForFile = ['window_',num2str(auxDataTest.readLength),'_withNor_MergeForwardAndReverse'];


clear fileName inputName outputName
for i=2:length(basicName)
  fileName{i} = ['data_',auxDataTest.basicDirNameForFile,'_normReadsf_',basicName{i}];
  outputName{i} = ['test_',auxDataTest.basicDirNameForFile,'_normReadsf_',basicName{i}];
end

if 1==2 % normReadsf
  userDir = getuserdir;
  % multiply the reads by 10000
  for i=2:length(basicName)
    
    % normalized Window90Norm_S1.mat
    dat = load([userDir,'/CS/BAC/12Samples/Solexa/data/',auxDataTest.basicDirNameForFile,'/Window90Norm_',basicName{i}]);
    
    datUnique = load([userDir,'/CS/BAC/12Samples/Solexa/data/window_90_noNor_MergeForwardAndReverse/data_window_90_noNor_MergeForwardAndReverse_',basicName{i}]);
    clear uniqueReads uniquereads_length
    uniqueReads = datUnique.uniqueReads;
    uniqueReads_length = ceil(dat.normReadsf*10^5);
    if size(uniqueReads_length,1)>size(uniqueReads_length,2)
      uniqueReads_length = uniqueReads_length';
    end
    length(find(uniqueReads_length==1))
    
    save([userDir,'/CS/BAC/12Samples/Solexa/data/',auxDataTest.basicDirNameForFile,'/',fileName{i}],'uniqueReads','uniqueReads_length')
  end
end

%%%%%%%%%%%%%%%%


% the same for normReads
%%%%% for normReadsf
clear
basicName = {'S1','S2','S3','S4','M1','M2','M3','M4','S7','S10','O7','O10'};

auxDataTest.datasetName = 'primers750_primer_tail';
auxDataTest.numBACtoConsider = 212040;
auxDataTest.dataSetPrimers750Flag = 1;

auxDataTest.readLength = 90;

auxDataTest.currNumProcessors = 20;
auxDataTest.basicDirNameForFile = ['window_',num2str(auxDataTest.readLength),'_withNor_MergeForwardAndReverse'];


clear fileName inputName outputName
for i=2:length(basicName)
  fileName{i} = ['data_',auxDataTest.basicDirNameForFile,'_normReads_',basicName{i}];
  outputName{i} = ['test_',auxDataTest.basicDirNameForFile,'_normReads_',basicName{i}];
end

if 1==2 % normReads
  userDir = getuserdir;
  % multiply the reads by 10000
  for i=2:length(basicName)
    
    % normalized Window90Norm_S1.mat
    dat = load([userDir,'/CS/BAC/12Samples/Solexa/data/',auxDataTest.basicDirNameForFile,'/Window90Norm_',basicName{i}]);
    
    datUnique = load([userDir,'/CS/BAC/12Samples/Solexa/data/window_90_noNor_MergeForwardAndReverse/data_window_90_noNor_MergeForwardAndReverse_',basicName{i}]);
    clear uniqueReads uniquereads_length
    uniqueReads = datUnique.uniqueReads;
    uniqueReads_length = ceil(dat.normReads*10^5);
    if size(uniqueReads_length,1)>size(uniqueReads_length,2)
      uniqueReads_length = uniqueReads_length';
    end
    
    length(find(uniqueReads_length==1))
    
    save([userDir,'/CS/BAC/12Samples/Solexa/data/',auxDataTest.basicDirNameForFile,'/',fileName{i}],'uniqueReads','uniqueReads_length')
  end
end


