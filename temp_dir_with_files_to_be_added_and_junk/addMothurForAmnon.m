function addMothurForAmnon(userDir,correctWeight,found,tmpMixtureName)
% mkdir ~/CS/BAC/mothurResFor6Region
% correctWeight = zeros(1,410849);correctWeight(1:10) = 0.1;found = zeros(1,410849);found(2:11) = 0.1;addMothurForAmnon('~',correctWeight,found,'j1')

%keyboard
%disp('note that the database used is the 410849')
%pause

userDirForFile = userDir;
%keyboard

auxData = struct;
auxData.basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/']; 
auxData.basicSeqKey = [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
outFileDirName = ['test_',tmpMixtureName];
auxData.tmpFileName = outFileDirName;

auxData.currDir = [userDirForFile,'/CS/BAC/tmpRuns/',outFileDirName];
auxData.saveName = [userDirForFile,'/CS/BAC/',outFileDirName];
auxData.userDir = userDir;

unix(['mkdir ',auxData.currDir]);

cd('~/')
pwd1 = pwd;
cd(auxData.currDir)
saveToFileName = [auxData.currDir,'/mothurRes_',tmpMixtureName];
CreateMothurDistNoam6Regions(saveToFileName,found,correctWeight,auxData);
cd(pwd1);
%keyboard
unix(['\mv -f ',auxData.currDir,'/mothurRes_',tmpMixtureName,'_mothurResL1L2.mat ',userDirForFile,...
      '/CS/BAC/mothurResFor6Region'])

unix(['rm -fr ',userDirForFile,'/CS/BAC/tmpRuns/',outFileDirName])
