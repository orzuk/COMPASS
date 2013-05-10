function addMothur(userDir,saveName,currSet,mixtureName,Nreads)
%keyboard
load(saveName)

tmpInd = store_kp{end};

nonEmpty = find(cellfun(@isempty,found)==0);
nonZer = find(abs(found{nonEmpty(end)})>10^-5);

correctWeight = zeros(1,410849);
correctWeight(currSet) = 0.01;

userDirForFile = userDir;
%keyboard

auxData = struct;
auxData.basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/']; 
auxData.basicSeqKey = [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
outFileDirName = ['test_',mixtureName,'_reads_readlen_100_noise_1_Nreads_1000000_correction_1_withCorrection_readlen_100'];
auxData.tmpFileName = outFileDirName;
auxData.currDir = [userDirForFile,'/CS/BAC/tmpRuns/',outFileDirName,'/',auxData.tmpFileName];
auxData.currDir = [userDirForFile,'/CS/BAC/tmpRuns/',outFileDirName,'/',auxData.tmpFileName];
auxData.saveName = [userDirForFile,'/CS/BAC/',outFileDirName,'/',auxData.tmpFileName];
auxData.userDir = userDir;

unix(['mkdir ',userDirForFile,'/CS/BAC/tmpRuns/',outFileDirName])
unix(['mkdir ',auxData.currDir])
unix(['mkdir ',userDirForFile,'/CS/BAC/',outFileDirName])

cd('~/')
pwd1 = pwd;
cd(auxData.currDir)
CreateMothurDistNoamL1L2(auxData.saveName,[],[],nonZer',abs(found{nonEmpty(end)}(nonZer))',correctWeight',auxData);
cd(pwd1);
unix(['\mv -f ',userDirForFile,'/CS/BAC/',outFileDirName,...
      '/test_',mixtureName,'_reads_readlen_100_noise_1_Nreads_',sprintf('%d',Nreads),'_correction_1_withCorrection_readlen_100_mothurResL1L2.mat ',userDirForFile,...
      '/CS/BAC/toyExample/EffectOfCloseBac/'])

unix(['rm -fr ',userDirForFile,'/CS/BAC/tmpRuns/',outFileDirName,';'...
      'rmdir ',userDirForFile,'/CS/BAC/tmpRuns/',outFileDirName,';'...
     'rmdir ',userDirForFile,'/CS/BAC/',outFileDirName])