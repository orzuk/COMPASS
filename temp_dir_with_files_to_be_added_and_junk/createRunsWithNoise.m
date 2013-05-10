% ff
clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
addpath([userDir,'/CS/mFilesBAC'])
outFileDirName = 'ff';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName])
unix(['mkdir ',userDir,'/CS/BAC/',outFileDirName])
prefix = outFileDirName;
brFlag = 1;

crGeneral(prefix,outFileDirName,20,[0],[10 100 100],[0 0.5],10,[inf],0,0,[50 100 400],'hour','week',brFlag);

name = 'tr1_51.txt'
fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name],'w');
fprintf(fid,'#!/bin/bash \n');
k = 1;
for i=[1:10:51]
  w = ['./listOfR_',prefix,'_',num2str(i),';sleep 20;'];
  fprintf(fid,'%s\n',w);
  if mod(k,10)==0
    fprintf(fid,'sleep 10m\n');
  end
  k = k+1;
end
fclose(fid);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name])

%%%%%%%%%%%%%%
clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
addpath([userDir,'/CS/mFilesBAC'])
outFileDirName = 'ff2';
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


crGeneral(prefix,outFileDirName,10,[0],[10 100 1000],[0 0.5],100,[inf],0,0,[50 100 400],'hour','week',brFlag,1);
crGeneral(prefix,outFileDirName,10,[0],[10 100 1000],[0 0.5],100,[inf],0,0,[50 100 400],'hour','week',brFlag,0);

name = 'tr1.txt'
fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name],'w');
fprintf(fid,'#!/bin/bash \n');
k = 1;
for i=[671:10:1791]
  w = ['./listOfR_',prefix,'_',num2str(i),';sleep 20;'];
  fprintf(fid,'%s\n',w);
  if mod(k,6)==0
    fprintf(fid,'sleep 30m\n');
  end
  k = k+1;
end
fclose(fid);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name])


crGeneral(prefix,outFileDirName,10,[0],[500],[0 0.5],100,[inf],0,0,[50 100 400],'hour','week',brFlag,1);
crGeneral(prefix,outFileDirName,10,[0],[500],[0 0.5],100,[inf],0,0,[50 100 400],'hour','week',brFlag,0);

name = 'tr2.txt'
fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name],'w');
fprintf(fid,'#!/bin/bash \n');
k = 1;
for i=[1:10:591]
  w = ['./listOfR_',prefix,'_',num2str(i),';sleep 20;'];
  fprintf(fid,'%s\n',w);
  if mod(k,6)==0
    fprintf(fid,'sleep 30m\n');
  end
  k = k+1;
end
fclose(fid);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name])

%%%%%%%%%%%%%%%%%%
% check the effect of averaging in repeatRandomGroups
clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
addpath([userDir,'/CS/mFilesBAC'])
outFileDirName = 'ffMeanInsteadOfMajority';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName])
unix(['mkdir ',userDir,'/CS/BAC/',outFileDirName])
prefix = outFileDirName;
brFlag = 1;


crGeneral(prefix,outFileDirName,10,[0],[10 100 1000],[0 0.5],100,[inf],0,0,[50 100 400],'hour','week',brFlag,1);
crGeneral(prefix,outFileDirName,10,[0],[10 100 1000],[0 0.5],100,[inf],0,0,[50 100 400],'hour','week',brFlag,0);


name = 'tr1.txt'
fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name],'w');
fprintf(fid,'#!/bin/bash \n');
k = 1;
for i=[1:10:691]
  w = ['./listOfR_',prefix,'_',num2str(i),';sleep 20;'];
  fprintf(fid,'%s\n',w);
  if mod(k,6)==0
    fprintf(fid,'sleep 30m\n');
  end
  k = k+1;
end
fclose(fid);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name])



%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% third
clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
addpath([userDir,'/CS/mFilesBAC'])
outFileDirName = 'ffThird';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/',outFileDirName])
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

crGeneralThird(prefix,outFileDirName,10,[0],[100],[0 0.5],20,[inf],0,0,[50 100 400],'hour','week',150000,20000,400,brFlag,1);
crGeneralThird(prefix,outFileDirName,10,[0],[100],[0 0.5],20,[inf],0,0,[50 100 400],'hour','week',150000,20000,400,brFlag,0);


name = 'tr1.txt'
fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name],'w');
fprintf(fid,'#!/bin/bash \n');
k = 1;
for i=[1:10:111]
  w = ['./listOfR_',prefix,'_',num2str(i),';sleep 20;'];
  fprintf(fid,'%s\n',w);
  if mod(k,6)==0
    fprintf(fid,'sleep 90m\n');
  end
  k = k+1;
end
fclose(fid);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name])

basicDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/ffThird/';
clear res
for readLength=[50 100 400]
  for np=[0 0.5]
    res(readLength,np*10+1,:) = extractRes(basicDir,'ffThird',0,100,np,Inf,0,0,readLength,410849);
  end
end



%%%%%%%%%%%%%%%%%%
clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
addpath([userDir,'/CS/mFilesBAC'])
outFileDirName = 'ffThirdTwiceLower20000';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/',outFileDirName])
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

crGeneralThird(prefix,outFileDirName,10,[0],[10 100 1000],[0 0.5],100,[inf],0,0,[50 100 400],'hour','week',150000,20000,400,brFlag,1);
crGeneralThird(prefix,outFileDirName,10,[0],[10 100],[0 0.5],100,[inf],0,0,[50 100 400],'hour','week',150000,20000,400,brFlag,0);

crGeneralThird(prefix,outFileDirName,10,[0],[10],[0 0.5],100,[inf],0,0,[700],'hour','week',150000,20000,400,brFlag,1);
crGeneralThird(prefix,outFileDirName,10,[0],[10],[0 0.5],100,[inf],0,0,[700],'hour','week',150000,20000,400,brFlag,0);

crGeneralThird(prefix,outFileDirName,10,[0],[100],[0 0.5],100,[inf],0,0,[700],'hour','week',150000,20000,400,brFlag,1);
crGeneralThird(prefix,outFileDirName,10,[0],[100],[0 0.5],100,[inf],0,0,[700],'hour','week',150000,20000,400,brFlag,0);

crGeneralThird([prefix,'_2'],outFileDirName,10,[0],[1000],[0 0.5],100,[inf],0,0,[50 100 400 700],'hour','week',150000,20000,400,brFlag,1);
crGeneralThird([prefix,'_2'],outFileDirName,10,[0],[1000],[0 0.5],100,[inf],0,0,[50 100 400 700],'hour','week',150000,20000,400,brFlag,0);

[junk,list50] = unix(['grep -l "readlen_50" ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/list*']);
[junk,list100] = unix(['grep -l "readlen_100" ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/list*']);
[junk,list400] = unix(['grep -l "readlen_400" ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/list*']);
[junk,list700] = unix(['grep -l "readlen_700" ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/list*']);


name = 'tr1_2.txt'
fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name],'w');
fprintf(fid,'#!/bin/bash \n');

listThem(fid,'120m',list50)
listThem(fid,'60m',list100)
listThem(fid,'40m',list400)
listThem(fid,'40m',list700)
fclose(fid);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name])



userDir = getuserdir;
%userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
basicDir = [userDir,'/CS/BAC/',outFileDirName,'/'];

outFileDirName = 'ffThirdTwiceLower20000_2';
clear res
nb = 1000;
for readLength=[50 100 400]
  for np=[0 0.5]
    res(readLength,np*10+1,:) = extractRes(basicDir,outFileDirName,0,nb,np,Inf,0,0,readLength,410849,1);
  end
end

%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%
clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
addpath([userDir,'/CS/mFilesBAC'])
outFileDirName = 'ffFourth';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/',outFileDirName])
prefix = outFileDirName;
brFlag = 1;

crGeneralThird(prefix,outFileDirName,10,[0],[10 100 1000],[0 0.5],100,[inf],0,0,[50 100 400],'hour','week',150000,20000,400,brFlag,1);
crGeneralThird(prefix,outFileDirName,10,[0],[10 100],[0 0.5],100,[inf],0,0,[50 100 400],'hour','week',150000,20000,400,brFlag,0);

[junk,list50] = unix(['grep -l "readlen_50" ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/list*']);
[junk,list100] = unix(['grep -l "readlen_100" ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/list*']);
[junk,list400] = unix(['grep -l "readlen_400" ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/list*']);
[junk,list700] = unix(['grep -l "readlen_700" ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/list*']);


name = 'tr1_2.txt'
fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name],'w');
fprintf(fid,'#!/bin/bash \n');

listThem(fid,'120m',list50)
listThem(fid,'60m',list100)
listThem(fid,'40m',list400)
listThem(fid,'40m',list700)
fclose(fid);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name])

%%%%%%%%%%%%%%%%%%%




% maybe randomize order
% change setting according to num_bac_in_mix

% check time according to num_bac_in_mix

% started 13:45

cd ~/CS/tmp; tar cf ffPriority.tar ffPriority;cd ~/CS/BAC/dataForSim/; tar cf ffPriority.tar ffPriority;cd ~/CS/BAC/tmpRuns/; tar cf ffPriority.tar ffPriority;

scp ~/CS/tmp/ffPriority.tar orzuk@tin.broadinstitute.org:/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp
scp ~/CS/BAC/dataForSim/ffPriority.tar  orzuk@tin.broadinstitute.org:/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/
scp ~/CS/BAC/tmpRuns/ffPriority.tar orzuk@tin.broadinstitute.org:/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/tmpRuns/

% in Br
cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp; tar xf ffPriority.tar;rm ffPriority.tar;cd /seq/orzuk2/compressed_sensing/ ...
    metagenomics/next_gen/CS/BAC/dataForSim/; tar xf ffPriority.tar;rm ffPriority.tar;cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/ ...
    BAC/tmpRuns/; tar xf ffPriority.tar;rm ffPriority.tar
