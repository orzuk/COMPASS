clear

% 10

basicDir = '~/CS/BAC/ff2/';
clear res
for readLength=[50,100,400]
  for np=[0 0.5]
    res(readLength,np*10+1,:) = extractRes(basicDir,'ff2',0,10,np,Inf,0,0,readLength,410849);
  end
end

basicDir = '~/CS/BAC/ffMeanInsteadOfMajority/';
clear res
for readLength=[50,100,400]
  for np=[0 0.5]
    res(readLength,np*10+1,:) = extractRes('~/CS/BAC/ffMeanInsteadOfMajority/','ffMeanInsteadOfMajority',0,10,np,Inf,0,0,readLength,410849);
  end
end


errorbar([50,100,400],res([50,100,400],1,1),res([50,100,400],1,2))
hold on
errorbar([50,100,400]+10,res([50,100,400],6,1),res([50,100,400],6,2),'color','r')

plot([50,100,400],res([[50,100,400],1]),'.')

% 25 75

% send 100 and 1000


%%%%%%%%%%%%%%%%%%55
dirName = 'ff3';
nameToTest = 'ff3_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_72_Nreads_Inf_Noise_0_Correction_0';


[w,addName] = nameAndAddName(dirName,nameToTest);



% run 
load(['~/tmp/structFor_',nameToTest]) 
%load(['~/tmp/',nameToTest,addName])
userdir = getuserdir;
auxData.saveName = [''];
auxData.saveName = [''];
auxData.currDir = [''];
auxData.basicSeqNameDir = [userdir,'/CS/BAC/datNoNonACGT/packed64/'];
auxData.basicSeqKey = [userdir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
auxData.reads_data_fileName = [userdir,'/tmp/',nameToTest,addName]
  
auxData.numProcessors = 7;
auxData.brFlag = 0;
  
save(['~/CS/tmp/run_',nameToTest],'auxData')
  
distributeBAC_generalOrSecond(['~/CS/tmp/run_',nameToTest]);

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
load ~/tmp/tmp_ff1
load ~/tmp/ffMeanInsteadOfMajority_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_2_Nreads_Inf_Noise_0_Correction_0_Nbacmix_10_Nread_Inf_Readlen_50_Npower_05_bacdistflag_0_NoReads
clear tmpInd1 tmpInd
k = 1;
part = 1:1000:length(store_kp{2});
part(end) = length(store_kp{2})+1;
for j=1:50
  A = randperm(length(store_kp{2}));
  for i=1:length(part)-1
    tmpInd1{k} = part(i):part(i+1)-1;
    tmpInd1{k} = store_kp{2}(A(tmpInd1{k}));
    k = k+1;
  end
end

matlabpool open local 7
parfor j=1:2500
    tmpI = tmpInd1(j);
    [junk,i1,i2] = intersect(tmpI{1},ind_bac_in_mix);
    if ~isempty(i1)
      [xx,sumRelevantReads]=solveForGroupOr(1,1,auxData,[],[],tmpI,auxData.basicSeqNameDir,auxData.basicSeqKey,auxData.readLength,uniqueReads,uniqueReads_length);
      if ~isempty(i1) && ~isempty(find(xx{1}(i1)<10^-3))
        j
        disp('found')
        save(['~/tmp/resMiss_',num2str(j)],'tmpI')
        %pause
      end
      xx{1}(i1)
    end
   
end
matlabpool close

  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
cd ~/tmp
load ffMeanInsteadOfMajority_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_Inf_Noise_0_Correction_0_Nbacmix_10_Nread_Inf_Readlen_50_Npower_05_bacdistflag_0_NoReads 
load structFor_ffMeanInsteadOfMajority_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_1_Nreads_Inf_Noise_0_Correction_0


clear tmpInd;

k = 1;
for j=1:100
  A = randperm(length(store_kp{2}));
  A = store_kp{2}(A(1:1000));
  A = setdiff(A,ind_bac_in_mix);
  for i=1:length(ind_bac_in_mix)
    tmpInd{k} = [ind_bac_in_mix(i);A'];
    k = k+1;
  end
end

matlabpool open local 7
parfor j=15:1000
    tmpI = tmpInd(j);
    [xx,sumRelevantReads]=solveForGroupOr(1,1,auxData,[],[],tmpI,auxData.basicSeqNameDir,auxData.basicSeqKey,auxData.readLength,uniqueReads,uniqueReads_length);
    
    [junk,i1,i2] = intersect(tmpI{1},ind_bac_in_mix);
    if ~isempty(find(xx{1}(i1)<10^-3))
      i
      disp('found')
      save(['~/tmp/resMiss_',num2str(j)],'tmpI')
      %pause
    end
    xx{1}(i1)
end
matlabpool close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% another way
clear tmpInd1 tmpInd
k = 1;
part = 1:1000:length(A);
part(end) = length(A)+1;
for j=1:50
  A = randperm(length(store_kp{2}));
  for i=1:length(part)-1
    tmpInd1{k} = part(i):part(i+1)-1;
    tmpInd1{k} = store_kp{2}(A(tmpInd1{k}));
    k = k+1;
  end
end

matlabpool open local 7
parfor j=1:2600
    tmpI = tmpInd1(j);
    [junk,i1,i2] = intersect(tmpI{1},ind_bac_in_mix);
    if ~isempty(i1)
      [xx,sumRelevantReads]=solveForGroupOr(1,1,auxData,[],[],tmpI,auxData.basicSeqNameDir,auxData.basicSeqKey,auxData.readLength,uniqueReads,uniqueReads_length);
      if ~isempty(i1) && ~isempty(find(xx{1}(i1)<10^-3))
        j
        disp('found')
        save(['~/tmp/resMiss_',num2str(j)],'tmpI')
        %pause
      end
      xx{1}(i1)
    end
   
end
matlabpool close

%%%%%%%%%%%%%%%%%%%%%%%55
clear
userDir = getuserdir;
addpath([userDir,'/CS/mFilesBAC'])
outFileDirName = 'ff1';
unix(['mkdir ',userDir,'/CS/tmp/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/dataForSim/',outFileDirName,...
      ';mkdir ',userDir,'/CS/BAC/tmpRuns/',outFileDirName])
unix(['mkdir ',userDir,'/CS/BAC/',outFileDirName])
prefix = outFileDirName;
brFlag = 0;

crGeneral(prefix,outFileDirName,7,[0],[10],[0],1,[inf],0,0,[400],'hour','week',brFlag,0);





%%%%%%%%%%%%%%%%%%%%%%%%%






basicDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/ff3/';
clear res
for readLength=[50]
  for np=[0.5]
    res(readLength,np*10+1,:) = extractRes(basicDir,'ff3',0,10,np,Inf,0,0,readLength,410849);
  end
end


clear
dirName = 'ff3';
nameToTest = 'ff3_bac_dist_flag_0_Nbac_in_mixture_10_npower_5_readlen_50_numIter_72_Nreads_Inf_Noise_0_Correction_0';

[w,addName] = nameAndAddName(dirName,nameToTest);


load(['~/tmp/structFor_',nameToTest]) 
load ~/tmp/dataForAnalysis_ff3
load(['~/tmp/',nameToTest,addName])

userdir = getuserdir;
auxData.saveName = [''];
auxData.saveName = [''];
auxData.currDir = [''];
auxData.basicSeqNameDir = [userdir,'/CS/BAC/datNoNonACGT/packed64/'];
auxData.basicSeqKey = [userdir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];



[uniqueReads,uniqueReads_length,auxData.fracRelevantReadsForInfinity]=createReadsForInfiniteNumberOr(ind_bac_in_mix,correctWeight,auxData.readLength,auxData.basicSeqNameDir,auxData.basicSeqKey);

clear tmpInd1 tmpInd
k = 1;
part = 1:1000:length(r.store_kp{2});
part(end) = length(r.store_kp{2})+1;
for j=1:50
  A = randperm(length(r.store_kp{2}));
  for i=1:length(part)-1
    tmpInd1{k} = part(i):part(i+1)-1;
    tmpInd1{k} = r.store_kp{2}(A(tmpInd1{k}));
    k = k+1;
  end
end

matlabpool open local 7
parfor j=1:size(tmpInd1,2)
    tmpI = tmpInd1(j);
    [junk,i1,i2] = intersect(tmpI{1},r.ind_bac_in_mix);
    if ~isempty(i1)
      [xx,sumRelevantReads]=solveForGroupOr(1,1,auxData,[],[],tmpI,auxData.basicSeqNameDir,auxData.basicSeqKey,auxData.readLength,uniqueReads,uniqueReads_length);
      if ~isempty(i1) && ~isempty(find(xx{1}(i1)<10^-3))
        j
        disp('found')
        rnd = round(rand(1)*1000);
        save(['~/tmp/resMiss_',num2str(rnd)],'tmpI')
        %pause
      end
      xx{1}(i1)
    end
   
end
matlabpool close


%%%%%%%%%%%
clear
basicDir = '~/CS/BAC/ffMeanInsteadOfMajority/';
clear res
for readLength=[50]
  for np=[0]
    res(readLength,np*10+1,:) = extractRes(basicDir,'ffMeanInsteadOfMajority',0,100,np,Inf,0,0,readLength,410849,1);
  end
end

%%%%%%%%%%%%
basicDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/ff4/';
clear res
for readLength=[50 100 400]
  for np=[0 0.5]
    res(readLength,np*10+1,:) = extractRes(basicDir,'ff4',0,10,np,Inf,0,0,readLength,410849);
  end
end
for readLength=[50 100 400]
  for np=[0 0.5]
    res(readLength,np*10+1,:) = extractRes(basicDir,'ff4',0,100,np,Inf,0,0,readLength,410849);
  end
end



%%%%%%%%%%%%%%%%%


%%%%%
%
basicSeqNameDir = ['~/CS/BAC/datNoNonACGT/packed64/'];
basicSeqKey= ['~/CS/BAC/datNoNonACGT/keyNoNonACGT'];

load(basicSeqKey,'len_uni') 


z = 371680%344025%142940%410009%19145;%276628;
[Header1,Sequence1] = loadSeqNames(z,basicSeqNameDir,basicSeqKey);
seq = int2nt(unpack_seqs(Sequence1{1},len_uni(z),64));

z_c = find(x_correct)

clear ham
for j=1:length(z_c)
  [Header1,Sequence_c] = loadSeqNames(z_c(j),basicSeqNameDir,basicSeqKey);
  seq_c = int2nt(unpack_seqs(Sequence_c{1},len_uni(z_c(j)),64));
  [cscore,algn]=swalign(seq,seq_c,'Alphabet','NT');
  ham(j) = length(find(algn=='|'));
end


z_c = find(x_correct)


%z = 344025%193289;
z_other = find(x_found>10^-3);
z_other = setdiff(z_other,z_c);

clear ham
for j=1:length(z_other)
  [Header1,Sequence_c] = loadSeqNames(z_other(j),basicSeqNameDir,basicSeqKey);
  seq_c = int2nt(unpack_seqs(Sequence_c{1},len_uni(z_other(j)),64));
  [cscore,algn]=swalign(seq,seq_c,'Alphabet','NT');
  ham(j) = length(find(algn=='|'));
end

for j=z_other(find(ham==1500))
  [Header1,Sequence_c] = loadSeqNames(j,basicSeqNameDir,basicSeqKey);
  seq_c = int2nt(unpack_seqs(Sequence_c{1},len_uni(j),64));
  [cscore,algn]=swalign(seq,seq_c,'Alphabet','NT');
  
  
end

%%%%%%%%%%%%%%%%5555
% matrix size as a function of readLength

clear
basicSeqNameDir = ['~/CS/BAC/datNoNonACGT/packed64/'];
basicSeqKey= ['~/CS/BAC/datNoNonACGT/keyNoNonACGT'];

load(basicSeqKey,'len_uni') 

numIter = 100;
list_of_readLength = [50 100 400 700 1000]
sz = zeros(length(list_of_readLength),numIter);
num_bac_in_mix = 10;
for j=1:length(list_of_readLength)
  readLength = list_of_readLength(j);
  for i=1:numIter
    tmpInd{1} = randperm(410891);
    tmpInd{1} = tmpInd{1}(1:num_bac_in_mix);
    [normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOr(readLength,tmpInd{1},basicSeqNameDir,basicSeqKey);
    sz(j,i) = size(normalizedBac,1);
    j,sz(j,i)
  end
end

plot(mean(sz,2),'-*')

%%%%%%%%%%%%%%%%%%55
normalizedBac = [1 0.5;...
                 0 0.5];
fracRelevantReads = [2;1]/3;
               
normalizedBac = [1 1;...
                 0 1];
fracRelevantReads = [2;1];
               

only AA:
normalizedBac = [1 1;...
                 0 1];
fracRelevantReads = [2 0]';

only CAA
fracRelevantReads = [1 1]';

0.5 CAA 0.5 AA
fracRelevantReads = [2 1]';
