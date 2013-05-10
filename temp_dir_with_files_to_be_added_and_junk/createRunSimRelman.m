clear

load('~/CS/BAC/toyExample/simRelman/seqs')

lenSeq = zeros(1,size(seqs,1));
for i=1:length(lenSeq)
  lenSeq(i) = length(seqs{i});
end

reads_uni = char(55*ones(size(seqs,1),max(lenSeq)));
for i=1:length(lenSeq)
  reads_uni(i,1:lenSeq(i)) = seqs{i};
end


save ~/CS/BAC/toyExample/simRelman/seqsRelman reads_uni lenSeq

auxData = struct;
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
auxData.userDir = userDir;
auxData.blastPath = [auxData.userDir,'/CS/BAC/12Samples/454/idx/blast-2.2.17/bin/'];
dataName = 'seqsRelman';
auxData.basicName = ['test_',dataName];
auxData.basicDir = [auxData.userDir,'/'];
auxData.databaseName = [auxData.userDir,'/CS/BAC/primers750_primer_tail/fastaFile/bac16s_primers750_primer_tail_full_without_ambiguous.fa'];
auxData.blastPath = [auxData.userDir,'/CS/BAC/12Samples/454/idx/blast-2.2.17/bin/'];
auxData.saveFileName = dataName;
auxData.dataDir = [auxData.userDir,'/CS/BAC/toyExample/simRelman/'];
auxData.dataName = dataName;
auxData.numProcessors = 200;
auxData.queueName = 'hour';

unix(['mkdir ',auxData.userDir,'/CS/BAC/alignToDatabase/tmpRuns/alignTmp_',auxData.dataName])
unix(['mkdir ',auxData.userDir,'/CS/BAC/alignToDatabase/dataForSim/alignTmp_',auxData.dataName])
unix(['mkdir ',auxData.basicDir,';mkdir ',auxData.basicDir,'/tmp'])

alignToDatabse(auxData)


%myFindSeqInDB_generic(auxData.basicName,auxData.basicDir,auxData.databaseName,auxData.blastPath,1,2,auxData.saveFileName,auxData.dataDir,auxData.dataName)


%%%%%%%%%55
clear 
load ~/CS/BAC/toyExample/simRelman/seqs
load ~/CS/BAC/toyExample/simRelman/alignTmp_seqsRelman_finalRes
load ~/CS/BAC/primers750_primer_tail/bac16s_primers750_primer_tail_full_without_ambiguous.mat

H = cell(1,size(freqs,2));
freqForMixture = H;
clear num*
for i=1:size(freqs,2)
  currInd = find(freqs(:,i));
  tmpBAC = IDX(currInd);
  
  [vals1 inds1 num_dups] = get_duplicates2(tmpBAC);
  if vals1(1)~=0
    disp('error');pause
  end
  
  currFreq = zeros(length(vals1)-1,1);
  for j=2:length(vals1)
    currFreq(j-1) = sum(freqs(currInd(inds1{j}),i));
  end
  
  
  freqForMixture{i} = currFreq;
  freqForMixture{i} = freqForMixture{i}./sum(freqForMixture{i}); 
  vals1(1) = [];
  
  H{i}.header = cell(1,length(vals1));
  for j=1:length(vals1)
    H{i}.header{j} = Header_750_tail{vals1(j)};
  end
end

save ~/CS/BAC/toyExample/simRelman/headerAndFrequencyRelman H freqForMixture

%%%%%%%%%%%%%%
clear
load ~/CS/BAC/toyExample/simRelman/headerAndFrequencyRelman

mixtureNumber = 5;
sim_bac = H{mixtureNumber}.header;
curr_freq = freqForMixture{mixtureNumber};

Header750 = load('~/CS/BAC/primers750_primer_tail/bac16s_primers750_primer_tail_full_without_ambiguous','Header_750_tail');

clear posIn750
for i=1:length(sim_bac)
  i
  for j=1:length(Header750.Header_750_tail)
    for k=1:length(Header750.Header_750_tail{j})
      if ~isempty(findstr(Header750.Header_750_tail{j}{k},sim_bac{i}{1}))
        posIn750{i} = j;
      end
    end
    
  end
end

a = cell2mat(posIn750)
length(a)
length(unique(a))

if length(a)~=length(unique(a))
  disp('problem')
end


ind_bac_in_mix_750 = cell2mat(posIn750);
N_750 = length(Header750.Header_750_tail);
correctWeight_750 = zeros(1,N_750);
correctWeight_750(ind_bac_in_mix_750) = freqForMixture{mixtureNumber};
correctWeight_750 = correctWeight_750./sum(correctWeight_750);
% del those that are lower than 10^-3
a = find(correctWeight_750>0 & correctWeight_750<10^-3)
correctWeight_750(a) = 0;
correctWeight_750 = correctWeight_750./sum(correctWeight_750);



save(['~/CS/BAC/toyExample/simRelman/data_simRelman_750_mixture_',num2str(mixtureNumber)],'correctWeight_750') 



%%%%%%%%%%%%
% 350
clear 


load ~/CS/BAC/toyExample/simRelman/headerAndFrequencyRelman

mixtureNumber = 5;
sim_bac = H{mixtureNumber}.header;
curr_freq = freqForMixture{mixtureNumber};
Header350 = load('~/CS/BAC/primers_startAs750_length350/bac16s_1to350_primers450_SequenceAndHeader');


clear posIn350
for i=1:length(sim_bac)
  i
  for j=1:length(Header350.header_uni1to350_primers450)
    for k=1:length(Header350.header_uni1to350_primers450{j})
      if ~isempty(findstr(Header350.header_uni1to350_primers450{j}{k},sim_bac{i}{1}))
        posIn350{i} = j;
      end
    end
    
  end
end

notFound = [];
for i=1:length(posIn350)
  if isempty(posIn350{i})
    notFound = [notFound,i];
  end
end

f = 1:length(posIn350);
f(notFound) = [];

a = cell2mat(posIn350(f))
length(a)
length(unique(a))

if length(a)~=length(unique(a))
  disp('problem')
end


ind_bac_in_mix_350 = cell2mat(posIn350(f));
N_350 = length(Header350.header_uni1to350_primers450);
correctWeight_350 = zeros(1,N_350);
correctWeight_350(ind_bac_in_mix_350) = freqForMixture{mixtureNumber}(f);
correctWeight_350 = correctWeight_350./sum(correctWeight_350);
% del those that are lower than 10^-3
a = find(correctWeight_350>0 & correctWeight_350<10^-3)
correctWeight_350(a) = 0;
correctWeight_350 = correctWeight_350./sum(correctWeight_350);


save(['~/CS/BAC/toyExample/simRelman/data_simRelman_350_mixture_',num2str(mixtureNumber)],'correctWeight_350') 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% create reads for runs - no noise 


% 750 create reads for runs - no noise 
clear
mixtureNumber = 5;
load(['~/CS/BAC/toyExample/simRelman/data_simRelman_750_mixture_',num2str(mixtureNumber)])
forBlastFlag = 0;
userDir = getuserdir;
basicSeqNameDir = [userDir,'/CS/BAC/primers750/datNoNonACGT/packed64/'];
basicSeqKey = [userDir,'/CS/BAC/primers750/datNoNonACGT/keyNoNonACGT_primers750'];

readLength = 100;
addNoiseFlag = 0;
Nreads = 10^6;

dataName = ['Relman_',num2str(mixtureNumber),'_reads750_readlen_',num2str(readLength),'_noise_',num2str(addNoiseFlag),'_Nreads_',num2str(Nreads)];
saveFile = [userDir,'/CS/BAC/toyExample/simRelman/reads/',dataName];

runSpecificToyExample(correctWeight_750,Nreads,basicSeqKey,basicSeqNameDir,addNoiseFlag,readLength,saveFile,forBlastFlag);

% 750 create reads for runs - with noise 
clear
mixtureNumber = 5;
load(['~/CS/BAC/toyExample/simRelman/data_simRelman_750_mixture_',num2str(mixtureNumber)])

forBlastFlag = 0;
userDir = getuserdir;
basicSeqNameDir = [userDir,'/CS/BAC/primers750/datNoNonACGT/packed64/'];
basicSeqKey = [userDir,'/CS/BAC/primers750/datNoNonACGT/keyNoNonACGT_primers750'];


readLength = 100;
addNoiseFlag = 1;
Nreads = 10^6;

ErrorStruct = []; ErrorStruct.error_model = 'exponential';
ErrorStruct.baseline_error = 0.005; 
ErrorStruct.final_error = 0.03; 
                    
p_one_nuc_error_mat = [ -1   0.3  0.22  0.18
                    0.5    -1  0.22   0.6
                    0.35  0.15    -1  0.22
                    0.15  0.55  0.56   -1]; % , 16, length(p));
ErrorStruct.p_one_nuc_error_mat = p_one_nuc_error_mat;

dataName = ['Relman_',num2str(mixtureNumber),'_reads750_readlen_',num2str(readLength),'_noise_',num2str(addNoiseFlag),'_Nreads_',num2str(Nreads)];
saveFile = [userDir,'/CS/BAC/toyExample/simRelman/reads/',dataName];

runSpecificToyExample(correctWeight_750,Nreads,basicSeqKey,basicSeqNameDir,addNoiseFlag,readLength,saveFile,forBlastFlag,ErrorStruct);




% 350 create reads for runs - no noise 

clear
mixtureNumber = 5;
load(['~/CS/BAC/toyExample/simRelman/data_simRelman_350_mixture_',num2str(mixtureNumber)])
userDir = getuserdir;
forBlastFlag = 1;
basicSeqNameDir = [userDir,'/CS/BAC/primers_startAs750_length350/datNoNonACGT/packed64/'];
basicSeqKey= [userDir,'/CS/BAC/primers_startAs750_length350/datNoNonACGT/keyNoNonACGT_primers_startAs750_length350'];

load(['~/CS/BAC/toyExample/simRelman/data_simRelman_350_mixture_',num2str(mixtureNumber)]) 
readLength = 350;
addNoiseFlag = 0;
Nreads = 10^4;

dataName = ['Relman_',num2str(mixtureNumber),'_reads350_readlen_',num2str(readLength),'_noise_',num2str(addNoiseFlag),'_Nreads_',num2str(Nreads)];
saveFile = [userDir,'/CS/BAC/toyExample/simRelman/reads/',dataName];

runSpecificToyExample(correctWeight_350,Nreads,basicSeqKey,basicSeqNameDir,addNoiseFlag,readLength,saveFile,forBlastFlag)

%%%%%%%%%%%%%%%%% % end create reads no noise


% 350 create reads with noise

clear
mixtureNumber = 5;
load(['~/CS/BAC/toyExample/simRelman/data_simRelman_350_mixture_',num2str(mixtureNumber)])
userDir = getuserdir;
forBlastFlag = 1;
basicSeqNameDir = [userDir,'/CS/BAC/primers_startAs750_length350/datNoNonACGT/packed64/'];
basicSeqKey= [userDir,'/CS/BAC/primers_startAs750_length350/datNoNonACGT/keyNoNonACGT_primers_startAs750_length350'];

load(['~/CS/BAC/toyExample/simRelman/data_simRelman_750_mixture_',num2str(mixtureNumber)]) 
readLength = 350;
addNoiseFlag = 1;
Nreads = 10^4;
ErrorStruct = []; ErrorStruct.error_model = 'exponential';
ErrorStruct.baseline_error = 0.005; 
ErrorStruct.final_error = 0.03; 
                    
p_one_nuc_error_mat = [ -1   0.3  0.22  0.18
                    0.5    -1  0.22   0.6
                    0.35  0.15    -1  0.22
                    0.15  0.55  0.56   -1]; % , 16, length(p));
ErrorStruct.p_one_nuc_error_mat = p_one_nuc_error_mat;

dataName = ['Relman_',num2str(mixtureNumber),'_reads350_readlen_',num2str(readLength),'_noise_',num2str(addNoiseFlag),'_Nreads_',num2str(Nreads)];
saveFile = [userDir,'/CS/BAC/toyExample/simRelman/reads/',dataName];

runSpecificToyExample(correctWeight_350,Nreads,basicSeqKey,basicSeqNameDir,addNoiseFlag,readLength,saveFile,forBlastFlag,ErrorStruct)

%%%%%%%%%%%%%%%%% % end create reads no noise




% run no noise case

% run 750 no noise case
clear
mixtureNumber = 5;
if  ~isempty(findstr(getenv('HOSTNAME'),'openu'))
  userDir = getuserdir;
else
  userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
end

readLength = 100;
addNoiseFlag = 0;
Nreads = 10^6;

dataName = ['Relman_',num2str(mixtureNumber),'_reads750_readlen_',num2str(readLength),'_noise_',num2str(addNoiseFlag),'_Nreads_',num2str(Nreads)];


auxDataTest = struct;
auxDataTest.datasetName = 'primers750_primer_tail';
auxDataTest.numBACtoConsider = 212040;
auxDataTest.dataSetPrimers750Flag = 1;
auxDataTest.readLength = 100;
auxDataTest.currNumProcessors = 20; % change for the with correction case
auxDataTest.correctMeasurementForReadErrorsFlag = 0; % run with no correction

basicDirNameForFile = dataName;
outputName = ['test_',dataName];
fileName = dataName;

testLikeRealData_genericChooseCorrection(outputName,['toyExample/simRelman/reads/',basicDirNameForFile],auxDataTest)

unix(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_',outputName,'/run_sol_',outputName,'_noCorrection_readlen_',num2str(auxDataTest.readLength)])

unix(['cp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/sol_',outputName,'/sol_',outputName,'_noCorrection_readlen_',num2str(auxDataTest.readLength),'/sol_',outputName,'_noCorrection_readlen_',num2str(auxDataTest.readLength),'.mat',' /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/toyExample/simRelman/results/'])


%%%%%%%%
% run 350 no noise
clear
mixtureNumber = 5;

if  ~isempty(findstr(getenv('HOSTNAME'),'openu'))
  userDir = getuserdir;
else
  userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
end

readLength = 350;
addNoiseFlag = 0;
Nreads = 10^4;
dataName = ['Relman_',num2str(mixtureNumber),'_reads350_readlen_',num2str(readLength),'_noise_',num2str(addNoiseFlag),'_Nreads_',num2str(Nreads)];



auxData = struct;
auxData.dataName = dataName;
auxData.userDir = userDir;
auxData.basicDir = [auxData.userDir,'/'];
auxData.basicName = ['test_',dataName];
auxData.databaseName = [auxData.userDir,'/CS/BAC/12Samples/454/idx/old_fasta_files_without_ambiguous/454idx.fa'];
auxData.blastPath = [auxData.userDir,'/CS/BAC/12Samples/454/idx/blast-2.2.17/bin/'];
auxData.dataDir = [auxData.userDir,'/CS/BAC/toyExample/simRelman/reads'];
N = load([auxData.dataDir,'/',dataName]);
auxData.numProcessors = max([2,ceil(size(N.reads_uni,1)/100)]);
auxData.queueName = 'hour';

unix(['mkdir ',auxData.userDir,'/CS/BAC/alignToDatabase/tmpRuns/alignTmp_',auxData.dataName])
unix(['mkdir ',auxData.userDir,'/CS/BAC/alignToDatabase/dataForSim/alignTmp_',auxData.dataName])

alignToDatabse(auxData)

%myFindSeqInDB_generic(auxData.basicName,auxData.basicDir,auxData.databaseName,auxData.blastPath,1,2,'junk',auxData.dataDir,auxData.dataName)



% end run 350


% end run no noise case

%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%

% run with noise


%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%

% run 350 with noise
clear
mixtureNumber = 5;
if  ~isempty(findstr(getenv('HOSTNAME'),'openu'))
  userDir = getuserdir;
else
  userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
end

readLength = 350;
addNoiseFlag = 1;
Nreads = 10^4;
dataName = ['Relman_',num2str(mixtureNumber),'_reads350_readlen_',num2str(readLength),'_noise_',num2str(addNoiseFlag),'_Nreads_',num2str(Nreads)];

auxData = struct;
auxData.dataName = dataName;
auxData.userDir = userDir;
auxData.basicDir = [auxData.userDir,'/'];
auxData.basicName = ['test_',dataName];
auxData.databaseName = [auxData.userDir,'/CS/BAC/12Samples/454/idx/old_fasta_files_without_ambiguous/454idx.fa'];
auxData.blastPath = [auxData.userDir,'/CS/BAC/12Samples/454/idx/blast-2.2.17/bin/'];
auxData.dataDir = [auxData.userDir,'/CS/BAC/toyExample/simRelman/reads'];
N = load([auxData.dataDir,'/',dataName]);
auxData.numProcessors = ceil(size(N.reads_uni,1)/100);
auxData.queueName = 'hour';

unix(['mkdir ',auxData.userDir,'/CS/BAC/alignToDatabase/tmpRuns/alignTmp_',auxData.dataName])
unix(['mkdir ',auxData.userDir,'/CS/BAC/alignToDatabase/dataForSim/alignTmp_',auxData.dataName])
%unix(['mkdir ',auxData.basicDir,';mkdir ',auxData.basicDir,'/tmp'])

alignToDatabse(auxData)




% run 750 with noise case - no correction
clear
mixtureNumber = 5;

if  ~isempty(findstr(getenv('HOSTNAME'),'openu'))
  userDir = getuserdir;
else
  userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
end
readLength = 100;
addNoiseFlag = 1;
Nreads = 10^6;

auxDataTest = struct;
auxDataTest.datasetName = 'primers750_primer_tail';
auxDataTest.numBACtoConsider = 212040;
auxDataTest.dataSetPrimers750Flag = 1;
auxDataTest.readLength = 100;
auxDataTest.currNumProcessors = 20; % change for the with correction case
auxDataTest.correctMeasurementForReadErrorsFlag = 0; % run with correction

dataName = ['Relman_',num2str(mixtureNumber),'_reads750_readlen_',num2str(readLength),'_noise_',num2str(addNoiseFlag),'_Nreads_',num2str(Nreads)];

basicDirNameForFile = [dataName,'_correction_',num2str(auxDataTest.correctMeasurementForReadErrorsFlag)]
outputName = ['test_',basicDirNameForFile];

testLikeRealData_genericChooseCorrection(outputName,['toyExample/simRelman/reads/',dataName],auxDataTest)

unix(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_',outputName,'/run_sol_',outputName,'_noCorrection_readlen_',num2str(auxDataTest.readLength)])

unix(['cp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/sol_',outputName,'/sol_',outputName,'_noCorrection_readlen_',num2str(auxDataTest.readLength),'/sol_',outputName,'_noCorrection_readlen_',num2str(auxDataTest.readLength),'.mat',' /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/toyExample/simRelman/results/'])



% run 750 with noise case - with correction
clear
mixtureNumber = 5;

if  ~isempty(findstr(getenv('HOSTNAME'),'openu'))
  userDir = getuserdir;
else
  userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
end

readLength = 100;
addNoiseFlag = 1;
Nreads = 10^6;

auxDataTest = struct;
auxDataTest.datasetName = 'primers750_primer_tail';
auxDataTest.numBACtoConsider = 212040;
auxDataTest.dataSetPrimers750Flag = 1;
auxDataTest.readLength = 100;
auxDataTest.currNumProcessors = 50; % change for the with correction case
auxDataTest.correctMeasurementForReadErrorsFlag = 1; % run with correction

dataName = ['Relman_',num2str(mixtureNumber),'_reads750_readlen_',num2str(readLength),'_noise_',num2str(addNoiseFlag),'_Nreads_',num2str(Nreads)];

basicDirNameForFile = [dataName,'_correction_',num2str(auxDataTest.correctMeasurementForReadErrorsFlag)]
outputName = ['test_',basicDirNameForFile];

testLikeRealData_genericChooseCorrection(outputName,['toyExample/simRelman/reads/',dataName],auxDataTest)

unix(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_',outputName,'/run_sol_',outputName,'_withCorrection_readlen_',num2str(auxDataTest.readLength)])

unix(['cp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/sol_',outputName,'/sol_',outputName,'_withCorrection_readlen_',num2str(auxDataTest.readLength),'/sol_',outputName,'_withCorrection_readlen_',num2str(auxDataTest.readLength),'.mat',' /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/toyExample/simRelman/results/'])


% more:
% run with noise
% high number of reads, low number of reads
% lower number of reads in 454



 

 



