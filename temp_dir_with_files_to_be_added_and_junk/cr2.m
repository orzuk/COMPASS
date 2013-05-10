clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
addpath([userDir,'/CS/mFilesBAC'])
outFileDirName = 'ff3';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName])
unix(['mkdir ',userDir,'/CS/BAC/',outFileDirName])
prefix = outFileDirName;
brFlag = 1;

if 1==2
  unix(['rm -fr ',userDir,'/CS/tmp/',outFileDirName,';',...
        ';rm -fr ',userDir,'/CS/BAC/dataForSim/',outFileDirName,';',...
        ';rm -fr ',userDir,'/CS/BAC/tmpRuns/',outFileDirName,';',...
        ';rm -fr ',userDir,'/CS/BAC/',outFileDirName,';',...
        ';rmdir ',userDir,'/CS/tmp/',outFileDirName,';',...
        ';rmdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,';',...
        ';rmdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName,';',...
        ';rmdir ',userDir,'/CS/BAC/',outFileDirName,';'])
end

crGeneral(prefix,outFileDirName,10,[0],[10],[0.5],100,[inf],0,0,[50],'hour','week',brFlag,1);
crGeneral(prefix,outFileDirName,10,[0],[10],[0.5],100,[inf],0,0,[50],'hour','week',brFlag,0);



name = 'tr1.txt'
fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name],'w');
fprintf(fid,'#!/bin/bash \n');
k = 1;
for i=[1:10:91]
  w = ['./listOfR_',prefix,'_',num2str(i),';sleep 20;'];
  fprintf(fid,'%s\n',w);
  if mod(k,6)==0
    fprintf(fid,'sleep 30m\n');
  end
  k = k+1;
end
fclose(fid);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name])

%%%%%%%%%%%%%%%%%%555
% run with repeats in move than 60000
clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
addpath([userDir,'/CS/mFilesBAC'])
outFileDirName = 'ff3';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName])
unix(['mkdir ',userDir,'/CS/BAC/',outFileDirName])
prefix = outFileDirName;
brFlag = 1;
repeatWhenLowerThanThisValue = 60000;

crGeneral(prefix,outFileDirName,10,[0],[10],[0.5],100,[inf],0,0,[50],'hour','week',repeatWhenLowerThanThisValue,brFlag,1);
crGeneral(prefix,outFileDirName,10,[0],[10],[0.5],100,[inf],0,0,[50],'hour','week',repeatWhenLowerThanThisValue,brFlag,0);

name = 'tr1.txt'
fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name],'w');
fprintf(fid,'#!/bin/bash \n');
k = 1;
for i=[1:10:91]
  w = ['./listOfR_',prefix,'_',num2str(i),';sleep 20;'];
  fprintf(fid,'%s\n',w);
  if mod(k,6)==0
    fprintf(fid,'sleep 30m\n');
  end
  k = k+1;
end
fclose(fid);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name])
