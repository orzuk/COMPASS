%%%%%%%%%%%
% test
% for Amnon
clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
outFileDirName = 'testing';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName])
unix(['mkdir ',userDir,'/CS/BAC/',outFileDirName])
prefix = 's8';
crOrSecond(prefix,outFileDirName,100,[0],[200],[0],1,[10^5]);

auxDataFileName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/testing/structFor_s8_bac_dist_flag_0_Nbac_in_mixture_200_npower_0_readlen_50_numIter_1_Nreads_100000';
profile -memory on
distributeBAC_generalOrSecond(auxDataFileName);;profile report;profile off;;profile report;profile off;
%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%
% test
% for Amnon
clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
outFileDirName = 't2';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName])
unix(['mkdir ',userDir,'/CS/BAC/',outFileDirName])
prefix = 's8';
crOrSecond(prefix,outFileDirName,100,[0],[10,50,100:100:500],[0],50,[10^5]);

fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','rt.txt'],'w');
fprintf(fid,'#!/bin/bash \n');
k = 1;
for i=[1:10:341]
  w = ['./listOfR_',prefix,'_',num2str(i),';sleep 30;'];
  fprintf(fid,'%s\n',w);
  if mod(k,5)==0
    fprintf(fid,'sleep 1h\n');
  end
  k = k+1;
end
fclose(fid);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','rt.txt'])


bac_dist_flag = 0;
Nbac_in_mixture = [10,50,100:100:500];
Npower = 0;
Nreads = 10^5;
numIter = 50;
[res,resCell]=loadResOr('s8','t2',bac_dist_flag,Nbac_in_mixture,Npower,numIter,Nreads);

README_input = str2mat('dimension 1: bac_dist_flag: 1 are distant, 2 are close ones',...
        'dimension 2: number of bacteria Nbac_in_mixture',...
        'dimension 3: distribution. 1 are power law 0.5 and 2 is uniform',...
        'dimension 4: number of reads',...
        'dimension 5: are the different realizations');
README_output = ['saves the results up to 10^-5. The correctWeight is given as is. Each cell entry has 2-3 lines regularly. The first cell is the correctWeight. ' ...
          ']The other lines correspond to solutions with less than 1500 bacteria']

Nbac_for_name = changeForName(Nbac_in_mixture);
bac_dist_flag_for_name = changeForName(bac_dist_flag);
Npower_for_name = changeForName(Npower*10);
Nreads_for_name = changeForName(Nreads);

basicName = ['res_',outFileDirName,'_bac_dist_flag',...
             '_Nbac_in_mixture',Nbac_for_name,...
             '_npower',Npower_for_name,...
             '_readlen',Nreads_for_name,...
             '_numIter_',num2str(numIter),...
             '_Nreads',Nreads_for_name];

save([userDir,'/CS/BAC/',outFileDirName,'/',basicName],'bac_dist_flag','Nbac_in_mixture','Npower','Nreads','numIter','resCell','README_input','README_output');
                 


%%%%%%%%%%%%%%%%%555
clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
outFileDirName = 't3';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName])
unix(['mkdir ',userDir,'/CS/BAC/',outFileDirName])
prefix = outFileDirName;
crOrSecond(prefix,outFileDirName,100,[0],[500],[0],1,[Inf]);

auxDataFileName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/t3/structFor_t3_bac_dist_flag_0_Nbac_in_mixture_500_npower_0_readlen_50_numIter_1_Nreads_Inf'

load /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/t3/structFor_t3_bac_dist_flag_0_Nbac_in_mixture_500_npower_0_readlen_50_numIter_1_Nreads_Inf

save /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/t3/structFor_t3_bac_dist_flag_0_Nbac_in_mixture_500_npower_0_readlen_50_numIter_1_Nreads_Inf_50CPU

auxDataFileName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/t3/structFor_t3_bac_dist_flag_0_Nbac_in_mixture_500_npower_0_readlen_50_numIter_1_Nreads_Inf_50CPU';


profile -memory on
distributeBAC_generalOrSecond(auxDataFileName);;profile report;profile off;profsave
%Sat Jul  9 03:47:51 EDT 2011
17:33
bkill -q hour -u orzuk 0;bkill -q week -u orzuk 0
bkill -q hour -u orzuk 0;bkill -q week -u orzuk 0

% compare 10 run on 50 CPUs and 100 CPUS

clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
outFileDirName = 't4';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName])
unix(['mkdir ',userDir,'/CS/BAC/',outFileDirName])
prefix = outFileDirName;
crOrSecond(prefix,outFileDirName,50,[0],[500],[0],1,[Inf]);

auxDataFileName = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/t3/structFor_t3_bac_dist_flag_0_Nbac_in_mixture_500_npower_0_readlen_50_numIter_1_Nreads_Inf'


% 
profile on; distributeBAC_generalOrSecond(auxDataFileName);profile report;profile off;profsave





%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare variable number of reads
clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
outFileDirName = 't2';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName])
unix(['mkdir ',userDir,'/CS/BAC/',outFileDirName])
prefix = 's8';
crOrSecond(prefix,outFileDirName,100,[0 1],[10,100,500],[0],50,[10^4 10^5 10^6 inf]);

fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','rt.txt'],'w');
fprintf(fid,'#!/bin/bash \n');
k = 1;
for i=[1:10:341]
  w = ['./listOfR_',prefix,'_',num2str(i),';sleep 30;'];
  fprintf(fid,'%s\n',w);
  if mod(k,5)==0
    fprintf(fid,'sleep 1h\n');
  end
  k = k+1;
end
fclose(fid);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/','rt.txt'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% test doing it all
clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
outFileDirName = 't4';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName])
unix(['mkdir ',userDir,'/CS/BAC/',outFileDirName])
prefix = 's8';
crOrSecond(prefix,outFileDirName,100,[1 0],[10,50,100 500],[0 0.5],10,[10^4 10^5 10^6 Inf]);

name = 'tr1_631.txt'
fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name],'w');
fprintf(fid,'#!/bin/bash \n');
k = 1;
for i=[1:10:631]
  w = ['./listOfR_',prefix,'_',num2str(i),';sleep 20;'];
  fprintf(fid,'%s\n',w);
  if mod(k,10)==0
    fprintf(fid,'sleep 10m\n');
  end
  k = k+1;
end
fclose(fid);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name])


clear
outFileDirName = 't4';
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
bac_dist_flag = [1 0]
Nbac_in_mixture = [10,50,100 500];
Npower = [0 0.5];
numIter = 10;
Nreads = [10^4 10^5 10^6 Inf];
[resCell]=loadResOr2('s8','t4',bac_dist_flag,Nbac_in_mixture,Npower,numIter,Nreads);

README_input = str2mat('dimension 1: bac_dist_flag: 1 are distant, 2 are close ones',...
        'dimension 2: number of bacteria Nbac_in_mixture',...
        'dimension 3: distribution. 1 are power law 0.5 and 2 is uniform',...
        'dimension 4: number of reads',...
        'dimension 5: are the different realizations');
README_output = ['saves the results up to 10^-5. The correctWeight is given as is. Each cell entry has 2-3 lines regularly. The first cell is the correctWeight. ' ...
          ']The other lines correspond to solutions with less than 1500 bacteria']

bac_dist_flag_for_name = changeForName(bac_dist_flag);
Nbac_for_name = changeForName(Nbac_in_mixture);
bac_dist_flag_for_name = changeForName(bac_dist_flag);
Npower_for_name = changeForName(Npower*10);
Nreads_for_name = changeForName(Nreads);

basicName = ['res_',outFileDirName,...
             '_bac_dist_flag',bac_dist_flag_for_name,...
             '_Nbac_in_mixture',Nbac_for_name,...
             '_npower',Npower_for_name,...
             '_readlen',Nreads_for_name,...
             '_numIter_',num2str(numIter),...
             '_Nreads',Nreads_for_name];

save([userDir,'/CS/BAC/',outFileDirName,'/',basicName],'bac_dist_flag','Nbac_in_mixture','Npower','Nreads','numIter','resCell','README_input','README_output');
   


%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% test doing it all
clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
outFileDirName = 't5';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName])
unix(['mkdir ',userDir,'/CS/BAC/',outFileDirName])
prefix = 's8';
crOrSecond(prefix,outFileDirName,100,[1 0],[10,50,100 500],[0 0.5],100,[10^4 10^5 10^6 Inf]);

name = 'tr1_6391.txt'
fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name],'w');
fprintf(fid,'#!/bin/bash \n');
k = 1;
for i=[641:10:6391]
  w = ['./listOfR_',prefix,'_',num2str(i),';sleep 20;'];
  fprintf(fid,'%s\n',w);
  if mod(k,10)==0
    fprintf(fid,'sleep 30m\n');
  end
  k = k+1;
end
fclose(fid);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name])



%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% test doing it all
clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
outFileDirName = 't6';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName])
unix(['mkdir ',userDir,'/CS/BAC/',outFileDirName])
prefix = 's8';
crOrSecond(prefix,outFileDirName,100,[1 0],[100 200 300 400 500],[0 0.5],10,[10^4 10^5 10^6 Inf]);

name = 'tr1_791.txt'
fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name],'w');
fprintf(fid,'#!/bin/bash \n');
k = 1;
for i=[1:10:791]
  w = ['./listOfR_',prefix,'_',num2str(i),';sleep 20;'];
  fprintf(fid,'%s\n',w);
  if mod(k,10)==0
    fprintf(fid,'sleep 10m\n');
  end
  k = k+1;
end
fclose(fid);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name])

saveToFile(prefix,outFileDirName,userDir,[1 0],[100 200 300 400 500],[0 0.5],10,[10^4 10^5 10^6 Inf])




%distributeBAC_generalOrSecond('/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/t4/structFor_s8_bac_dist_flag_0_Nbac_in_mixture_100_npower_0_readlen_50_numIter_1_Nreads_10000')
 

%%%%%%%%%
% test run with smaller number of CPU
clear
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
outFileDirName = 'test150_to_450_lowCPU';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName])
unix(['mkdir ',userDir,'/CS/BAC/',outFileDirName])
prefix = outFileDirName;
crOrSecond(prefix,outFileDirName,10,[1 0],[100:100:500],[0 0.5],10,[10^5 10^6]);

name = 'tr1_391.txt'
fid = fopen([userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name],'w');
fprintf(fid,'#!/bin/bash \n');
k = 1;
for i=[1:10:391]
  w = ['./listOfR_',prefix,'_',num2str(i),';sleep 20;'];
  fprintf(fid,'%s\n',w);
  if mod(k,10)==0
    fprintf(fid,'sleep 30m\n');
  end
  k = k+1;
end
fclose(fid);
unix(['chmod 0700 ',userDir,'/CS/BAC/dataForSim/',outFileDirName,'/',name])

