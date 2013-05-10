%%%%%%%%%%%%%%
% prepare the 90 reads

clear
basicName = {'S1','S2','S3','S4','M1','M2','M3','M4','S7','S10','O7','O10'};

auxDataTest = struct;

auxData = struct;
auxData.readLengthPost = 90;
auxData.readLengthPre = 100;
auxData.readLength = auxData.readLengthPost;
%auxData.inputDirName = '~/CS/BAC/12Samples/Solexa/data/opt3_th50_mar25_2012';
auxData.uniqueReadsDirName = '~/CS/BAC/12Samples/Solexa/data/opt3_MergeForwardReverse';
auxData.POS_dirName = '~/CS/BAC/findAlignment_MergeForwardAndReverse'; % new version
auxData.POS_basicFileName = 'alignTmp_MergeForwardAndReverse_100';
auxData.saveDir = ['~/CS/BAC/12Samples/Solexa/data/window_',num2str(auxData.readLengthPost),'_noNor_MergeForwardAndReverse'];

for i=2:12
  auxData.POS_fileName = [auxData.POS_basicFileName,'_',basicName{i}];
  auxData.uniqueReadsFilesName = ['data_NO_OPT3_opt3_MergeForwardReverse_',basicName{i}];
  auxData.fileName = ['data_window_',num2str(auxData.readLengthPost),'_noNor_MergeForwardAndReverse_',basicName{i}];

  % new version for MergeForwardAndReverse
  prepareWindow_12Samples_createReads_MergeForwardAndReverse(auxData)
  
end

%%%%%%%%%%%%%%%%
% end create reads


