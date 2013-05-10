%%%%%%%%%%%
% create
clear


cr('s5',[0],[500],[0 0.5],100,[10^4 10^5 Inf])

userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
fid = fopen([userDir,'/CS/BAC/dataForSim/set/sec/rt.txt'],'w');
fprintf(fid,'#!/bin/bash \n');
k = 1;
for i=[1:10:591]
  w = ['./listOfR_s5_',num2str(i),';sleep 30;'];
  fprintf(fid,'%s\n',w);
  if mod(k,5)==0
    fprintf(fid,'sleep 3h\n');
  end
  k = k+1;
end
fclose(fid);



%rerun s1
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
basicPrefix = 's5';
basicDir = '/CS/BAC/dataForSim/set/sec/';
basicRerunName = [basicPrefix,'_rerun1'];
load([userDir,basicDir,basicPrefix,'_list'])
rerun(totalListOfNames,totalListOfRuns,basicDir,basicPrefix,basicRerunName,userDir);


% del *.e files and run
unix(['cd ',userDir,basicDir,';rm ',basicPrefix,'*.e'])
unix(['cd ',userDir,basicDir,';chmod 0700 ',userDir,basicDir,basicRerunName])


% check the reads are the same
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
basicSeqNameDir = [userDir,'/CS/BAC/dat450000/'];
basicSeqKey =  [userDir,'/CS/BAC/dat450000/key450000'];

loadRes('s5',0,10,0.5,100,[Inf])
loadRes('s2',1,10,0,50,[10^5])


% 
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
basicFileName = [userDir,'/CS/BAC/sec/s5_bac_dist_flag_0_Nbac_in_mixture_500_npower_0'];
correctDataBasicFileName = [userDir,'/CS/BAC/dataForSim/set/sec/s5_bac_dist_flag_0_Nbac_in_mixture_500_npower_0_readlen_50_numIter_',];

unix(['tar cf ',userDir,'/rFiles.tar ',basicFileName,'*&']);
unix(['tar cf ',userDir,'/cFiles.tar ',correctDataBasicFileName,'*&']);
unix(['scp ',userDir,'/rFiles.tar shental@sol.cslab.openu.ac.il:~/CS/BAC/sec/'])
unix(['scp ',userDir,'/cFiles.tar shental@sol.cslab.openu.ac.il:~/CS/BAC/dataForSim/set/sec/'])


