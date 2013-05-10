%%%%%%%%%%%%%%%%%%555
% run with repeats in move than 60000
clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
addpath([userDir,'/CS/mFilesBAC'])
outFileDirName = 'ff4';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName])
unix(['mkdir ',userDir,'/CS/BAC/',outFileDirName])
prefix = outFileDirName;
brFlag = 1;
repeatWhenLowerThanThisValue = 60000;

crGeneral(prefix,outFileDirName,10,[0],[10 100 1000],[0 0.5],20,[inf],0,0,[50 100 400],'hour','week',repeatWhenLowerThanThisValue,brFlag,1);
crGeneral(prefix,outFileDirName,10,[0],[10 100 1000],[0 0.5],20,[inf],0,0,[50 100 400],'hour','week',repeatWhenLowerThanThisValue,brFlag,0);



name = 'tr1.txt'
fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name],'w');
fprintf(fid,'#!/bin/bash \n');
k = 1;
for i=[11:10:111]
  w = ['./listOfR_',prefix,'_',num2str(i),';sleep 20;'];
  fprintf(fid,'%s\n',w);
  if mod(k,4)==0
    fprintf(fid,'sleep 120m\n');
  end
  k = k+1;
end
fclose(fid);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name])
