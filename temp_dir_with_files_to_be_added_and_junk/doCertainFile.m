function doCertainFile(outFileDirName,name)

userDir = getuserdir;
name_L1Data = [userDir,'/CS/BAC/',outFileDirName,'/',name];
load(name_L1Data)
  
fileNameL1Res = [userDir,'/CS/BAC/',outFileDirName,'/',name];
a = findstr(fileNameL1Res,'_L1data');
fileNameL1Res(a:end) = [];
fileNameL1Res = [fileNameL1Res,'_L1Res']

addL1toL2_res(uniqueReads,uniqueReads_length,tmpInd,fileNameL1Res,auxData)
