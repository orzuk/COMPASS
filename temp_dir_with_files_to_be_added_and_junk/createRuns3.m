%%%%%%%%%%%
% create
clear

cr('s1',0,10,0.5,50,[Inf])
cr('s2',1,10,0,50,[10^5])




w = [];
for i=[11:10:41]
  w = [w,'./listOfR_s1_',num2str(i),';sleep 30;'];
end

w = [];
for i=[1:10:41]
  w = [w,'./listOfR_s2_',num2str(i),';sleep 30;'];
end

%rerun s1
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
basicPrefix = 's1';
basicDir = '/CS/BAC/dataForSim/set/sec/';
basicRerunName = [basicPrefix,'_rerun1'];
load([userDir,basicDir,basicPrefix,'_list'])
rerun(totalListOfNames,totalListOfRuns,basicDir,basicPrefix,basicRerunName,userDir);


% del *.e files and run
unix(['cd ',userDir,basicDir,';rm ',basicPrefix,'*.e'])x
unix(['cd ',userDir,basicDir,';chmod 0700 ',userDir,basicDir,basicRerunName])

%rerun s2
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
basicPrefix = 's2';
basicDir = '/CS/BAC/dataForSim/set/sec/';
basicRerunName = [basicPrefix,'_rerun1'];
load([userDir,basicDir,basicPrefix,'_list'])
rerun(totalListOfNames,totalListOfRuns,basicDir,basicPrefix,basicRerunName,userDir);


% del *.e files and run
unix(['cd ',userDir,basicDir,';rm ',basicPrefix,'*.e'])
unix(['cd ',userDir,basicDir,';chmod 0700 ',userDir,basicDir,basicRerunName])

%%%%%%%%
% move them
unix(['scp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/sec/s1* shental@sol.cslab.openu.ac.il:~/CS/BAC/sec/'])
unix(['scp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/sec/s2* shental@sol.cslab.openu.ac.il:~/CS/BAC/sec/'])

%%%%%%%%%%
% collrect the results

loadRes('s1',0,10,0.5,50,[Inf])
loadRes('s2',1,10,0,50,[10^5])

%%%%%%
% compare 
basicSeqKey = '/homes/csfaculty/shental/CS/BAC/dat450000/key450000';
basicSeqNameDir = '/homes/csfaculty/shental/CS/BAC/dat450000/';

resFile = '~/CS/BAC/sec/s1_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_Inf';
correctFile = '~/CS/BAC/dataForSim/set/sec/s1_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_Inf_Nbacmix_10_Nread_Inf_Readlen_50_Npower_05_bacdistflag_0_NoReads.mat';
dist = 51;
freqset=compareResults(dist,resFile,correctFile,basicSeqKey,basicSeqNameDir)

load(resFile)
load(correctFile)


plotRes(resFile,correctFile)

resFile = '~/CS/BAC/sec/s2_bac_dist_flag_1_Nbac_in_mixture_10_npower_0_readlen_50_numIter_9_Nreads_100000';
correctFile = '~/CS/BAC/dataForSim/set/sec/s2_bac_dist_flag_1_Nbac_in_mixture_10_npower_0_readlen_50_numIter_9_Nreads_100000_Nbacmix_10_Nread_100000_Readlen_50_Npower_0_bacdistflag_1_NoReads';
dist = 51;
freqset=compareResults(dist,resFile,correctFile,basicSeqKey,basicSeqNameDir)


plotRes('~/CS/BAC/sec/s2_bac_dist_flag_1_Nbac_in_mixture_10_npower_0_readlen_50_numIter_15_Nreads_100000','~/CS/BAC/dataForSim/set/sec/s2_bac_dist_flag_1_Nbac_in_mixture_10_npower_0_readlen_50_numIter_15_Nreads_100000_Nbacmix_10_Nread_100000_Readlen_50_Npower_0_bacdistflag_1_NoReads')


load ~/CS/BAC/sec/s2_bac_dist_flag_1_Nbac_in_mixture_10_npower_0_readlen_50_numIter_15_Nreads_100000
load ~/CS/BAC/dataForSim/set/sec/s2_bac_dist_flag_1_Nbac_in_mixture_10_npower_0_readlen_50_numIter_15_Nreads_100000_Nbacmix_10_Nread_100000_Readlen_50_Npower_0_bacdistflag_1_NoReads


basicSeqKey = '/homes/csfaculty/shental/CS/BAC/dat450000/key450000';
basicSeqNameDir = '/homes/csfaculty/shental/CS/BAC/dat450000/';

[normalizedBac values] = prepareGroupOf1000DistributedSequenceFiles2(50,store_kp{3},basicSeqNameDir,basicSeqKey);

[uniqueReads_Inf,uniqueReads_length_Inf,auxData.fracRelevantReadsForInfinity]=createReadsForInfiniteNumber(ind_bac_in_mix,correctWeight,50,basicSeqNameDir,basicSeqKey);

dataIn = struct;
dataIn.fracRelevantReadsForInfinity = auxData.fracRelevantReadsForInfinity;

[fracRelevantReads,sumRelevantReads] = currReads(uniqueReads_Inf,uniqueReads_length_Inf,values,1,dataIn);

[x]=runOneGroupOf1000ForCompilation(normalizedBac,fracRelevantReads,1);

s = zeros(size(correctWeight));
s(store_kp{3}) = x;

plot(s,correctWeight,'.')


% example results
for i=1:50
  unix(['cp ~/CS/BAC/sec/s2_bac_dist_flag_1_Nbac_in_mixture_10_npower_0_readlen_50_numIter_',num2str(i),'_Nreads_100000.mat ',['~/CS/BAC/' ...
                      'dataForSim/set/sec/s2_bac_dist_flag_1_Nbac_in_mixture_10_npower_0_readlen_50_numIter_'],num2str(i), ...
        '_Nreads_100000_Nbacmix_10_Nread_100000_Readlen_50_Npower_0_bacdistflag_1_NoReads.mat ',' ~/CS/BAC/files/'])
end


for i=1:50
  unix(['cp ~/CS/BAC/sec/s1_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_',num2str(i),'_Nreads_Inf.mat ',['~/CS/BAC/' ...
                      'dataForSim/set/sec/s1_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_'],num2str(i), ...
        '_Nreads_Inf_Nbacmix_10_Nread_Inf_Readlen_50_Npower_05_bacdistflag_0_NoReads.mat',' ~/CS/BAC/files/'])
end

