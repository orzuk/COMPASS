%%%%%%%%%%%
% create
clear

userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';

outFileDirName = 'testOr';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName])
crOr('s6',outFileDirName,100,[0],[500],[0],2,[Inf])



% test how to run
outFileDirName = 'testOrSmallNumCPU';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName,';mkdir ',userDir,'/CS/BAC/',outFileDirName])
crOr('s6',outFileDirName,11,[0],[500],[0],10,[Inf])

start at 11:00
outFileDirName = 'testOrLargeNumCPU';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName,';mkdir ',userDir,'/CS/BAC/',outFileDirName])
crOr('s6',outFileDirName,110,[0],[500],[0],1,[Inf])


fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','rt.txt'],'w');
fprintf(fid,'#!/bin/bash \n');
k = 1;
for i=[1]
  w = ['./listOfR_s5_',num2str(i),';sleep 30;'];
  fprintf(fid,'%s\n',w);
  if mod(k,5)==0
    fprintf(fid,'sleep 3h\n');
  end
  k = k+1;
end
fclose(fid);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','rt.txt'])



loadResOr('s6','testOrLargeNumCPU',[0],[500],[0],1,[Inf])
loadResOr('s6','testOr',[0],[500],[0],2,[Inf])


% del *.e files and run
unix(['cd ',userDir,basicDir,';rm ',basicPrefix,'*.e'])
unix(['cd ',userDir,basicDir,';chmod 0700 ',userDir,basicDir,basicRerunName])


% check the reads are the same
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
basicSeqNameDir = [userDir,'/CS/BAC/dat450000/'];
basicSeqKey =  [userDir,'/CS/BAC/dat450000/key450000'];



% 
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
basicFileName = [userDir,'/CS/BAC/sec/s5_bac_dist_flag_0_Nbac_in_mixture_500_npower_0'];
correctDataBasicFileName = [userDir,'/CS/BAC/dataForSim/set/sec/s5_bac_dist_flag_0_Nbac_in_mixture_500_npower_0_readlen_50_numIter_',];

unix(['tar cf ',userDir,'/rFiles.tar ',basicFileName,'*&']);
unix(['tar cf ',userDir,'/cFiles.tar ',correctDataBasicFileName,'*&']);
unix(['scp ',userDir,'/rFiles.tar shental@sol.cslab.openu.ac.il:~/CS/BAC/sec/'])
unix(['scp ',userDir,'/cFiles.tar shental@sol.cslab.openu.ac.il:~/CS/BAC/dataForSim/set/sec/'])


