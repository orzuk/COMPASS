
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%

% prepare the reads
clear
basicName = {'S1','S2','S3','S4','M1','M2','M3','M4','S7','S10','O7','O10'};

auxDataTest = struct;

auxData = struct;
auxData.readLengthPost = 90;
auxData.readLengthPre = 100;
auxData.readLength = auxData.readLengthPost;
%auxData.inputDirName = '~/CS/BAC/12Samples/Solexa/data/opt3_th50_mar25_2012';
auxData.uniqueReadsDirName = '~/CS/BAC/12Samples/Solexa/data/opt3_th50_mar25_2012';
auxData.POS_dirName = '~/CS/BAC/findAlignment_May12'; % new version
auxData.POS_fileName = 'alignTmp_May12_';
auxData.saveDir = ['~/CS/BAC/12Samples/Solexa/data/window_',num2str(auxData.readLengthPost),'_noNor'];

% sample related issues
% uses the [opt3_th50_mar25_2012/illumina_reads100_S1_unifreq_corr_opt3_th50]

% freq_read is the original read frequency before normalizing
% uni_read_seq are unique reads which are alos aligned - these are only forward and reverse compelement of reverse strand
% POS = is the alignment position. For reads originating from reverse strand - the POS is negative 

for i=1:12
  auxData.uniqueReadsFilesName = ['illumina_reads100_',basicName{i},'_unifreq_corr_opt3_th50'];
  auxData.fileName = ['data_window_',num2str(auxData.readLengthPost),'_noNor_',basicName{i}];
  % old version
  %prepareWindow_12Samples_createReads(auxData) % uses Amit's alignment data
  
  % new version with more data
  prepareWindow_12Samples_createReads_newVersion(auxData)
  
end
%%%%%%%%%%%%%%%%
% end create reads

%%%%%%%%%%%%%%%%%%%%
% find the position of the reads using "createRes_findPositionReads_May12"
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%
% 

basicDirNameForFile = ['window_',num2str(auxData.readLengthPost),'_noNor'];
% copy to br
PWD = pwd;
clear w
w = ['cd ~/CS/BAC/12Samples/Solexa/data/',basicDirNameForFile,';scp ']
for i=1:length(basicName)
  fileName = ['data_window_',num2str(auxData.readLengthPost),'_noNor_',basicName{i}];
  w = [w,' ',fileName,'.mat '];
end
w = [w,' orzuk@tin.broadinstitute.org:/seq/orzuk2/compressed_sensing/metagenomics/next_gen','/CS/BAC/12Samples/Solexa/data/',basicDirNameForFile];

cd(PWD);


