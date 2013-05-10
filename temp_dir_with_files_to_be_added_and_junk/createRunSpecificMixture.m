clear

fid=fopen('~/CS/BAC/toyExample/sim1/list_of_bacteria_for_toy_simulation');
k = 1;
tline = 1;
clear sim_bac
while ischar(tline) | tline~=-1
  if ischar(tline) && ~isempty(tline)
    sim_bac{k} = tline;
    k = k+1;
  end
  tline = fgetl(fid);
end
fclose(fid);
      
% delete 9 17 20 which did not appear in 750 (probably due to NaN)
sim_bac([9 17 20]) = []

% find the reads in the header file - 750 database
Header750 = load('~/CS/BAC/primers750_primer_tail/bac16s_primers750_primer_tail_full_without_ambiguous','Header_750_tail');

clear posIn750
for i=1:length(sim_bac)
  i
  for j=1:length(Header750.Header_750_tail)
    for k=1:length(Header750.Header_750_tail{j})
      if ~isempty(findstr(Header750.Header_750_tail{j}{k},sim_bac{i}))
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

ind_bac_in_mix_750 = unique(a);

N_750 = length(Header750.Header_750_tail);
correctWeight_750 = zeros(1,N_750);
correctWeight_750(ind_bac_in_mix_750) = 1/length(ind_bac_in_mix_750);

save ~/CS/BAC/toyExample/sim1/data_sim1_750 correctWeight_750



%%%%%%%%%%%%
% 350
Header350 = load('~/CS/BAC/primers_startAs750_length350/bac16s_1to350_primers450_SequenceAndHeader');
clear posIn350
for i=1:length(sim_bac)
  i
  for j=1:length(Header350.header_uni1to350_primers450)
    for k=1:length(Header350.header_uni1to350_primers450{j})
      if ~isempty(findstr(Header350.header_uni1to350_primers450{j}{k},sim_bac{i}))
        posIn350{i} = j;
      end
    end
    
  end
end

a = cell2mat(posIn350)
length(a)
length(unique(a))

ind_bac_in_mix_350 = unique(a);

N_350 = length(Header350.header_uni1to350_primers450);
correctWeight_350 = zeros(1,N_350);
correctWeight_350(ind_bac_in_mix_350) = 1/length(ind_bac_in_mix_350);

save ~/CS/BAC/toyExample/sim1/data_sim1_350 correctWeight_350

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create reads for runs - no noise 


% 750 create reads for runs - no noise 
clear
forBlastFlag = 0;
userDir = getuserdir;
basicSeqNameDir = [userDir,'/CS/BAC/primers750/datNoNonACGT/packed64/'];
basicSeqKey = [userDir,'/CS/BAC/primers750/datNoNonACGT/keyNoNonACGT_primers750'];

load ~/CS/BAC/toyExample/sim1/data_sim1_750
readLength = 100;
addNoiseFlag = 0;
Nreads = 10^6;

dataName = ['reads750_readlen_',num2str(readLength),'_noise_',num2str(addNoiseFlag),'_Nreads_',num2str(Nreads)];
saveFile = [userDir,'/CS/BAC/toyExample/sim1/reads/',dataName];

runSpecificToyExample(correctWeight_750,Nreads,basicSeqKey,basicSeqNameDir,addNoiseFlag,readLength,saveFile,forBlastFlag);

% 750 create reads for runs - with noise 
clear
forBlastFlag = 0;
userDir = getuserdir;
basicSeqNameDir = [userDir,'/CS/BAC/primers750/datNoNonACGT/packed64/'];
basicSeqKey = [userDir,'/CS/BAC/primers750/datNoNonACGT/keyNoNonACGT_primers750'];

load ~/CS/BAC/toyExample/sim1/data_sim1_750 
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

dataName = ['reads750_readlen_',num2str(readLength),'_noise_',num2str(addNoiseFlag),'_Nreads_',num2str(Nreads)];
saveFile = [userDir,'/CS/BAC/toyExample/sim1/reads/',dataName];

runSpecificToyExample(correctWeight_750,Nreads,basicSeqKey,basicSeqNameDir,addNoiseFlag,readLength,saveFile,forBlastFlag,ErrorStruct);




% 350 create reads for runs - no noise 

clear
userDir = getuserdir;
forBlastFlag = 1;
basicSeqNameDir = [userDir,'/CS/BAC/primers_startAs750_length350/datNoNonACGT/packed64/'];
basicSeqKey= [userDir,'/CS/BAC/primers_startAs750_length350/datNoNonACGT/keyNoNonACGT_primers_startAs750_length350'];

load ~/CS/BAC/toyExample/sim1/data_sim1_350
readLength = 350;
addNoiseFlag = 0;
Nreads = 10^4;

dataName = ['reads350_readlen_',num2str(readLength),'_noise_',num2str(addNoiseFlag),'_Nreads_',num2str(Nreads)];
saveFile = [userDir,'/CS/BAC/toyExample/sim1/reads/',dataName];

runSpecificToyExample(correctWeight_350,Nreads,basicSeqKey,basicSeqNameDir,addNoiseFlag,readLength,saveFile,forBlastFlag)

%%%%%%%%%%%%%%%%% % end create reads no noise


% 350 create reads with noise

clear
userDir = getuserdir;
forBlastFlag = 1;
basicSeqNameDir = [userDir,'/CS/BAC/primers_startAs750_length350/datNoNonACGT/packed64/'];
basicSeqKey= [userDir,'/CS/BAC/primers_startAs750_length350/datNoNonACGT/keyNoNonACGT_primers_startAs750_length350'];

load ~/CS/BAC/toyExample/sim1/data_sim1_350
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

dataName = ['reads350_readlen_',num2str(readLength),'_noise_',num2str(addNoiseFlag),'_Nreads_',num2str(Nreads)];
saveFile = [userDir,'/CS/BAC/toyExample/sim1/reads/',dataName];

runSpecificToyExample(correctWeight_350,Nreads,basicSeqKey,basicSeqNameDir,addNoiseFlag,readLength,saveFile,forBlastFlag,ErrorStruct)

%%%%%%%%%%%%%%%%% % end create reads no noise




% run no noise case

% run 750 no noise case
clear
if  ~isempty(findstr(getenv('HOSTNAME'),'openu'))
  userDir = getuserdir;
else
  userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
end

readLength = 100;
addNoiseFlag = 0;
Nreads = 10^6;

dataName = ['reads750_readlen_',num2str(readLength),'_noise_',num2str(addNoiseFlag),'_Nreads_',num2str(Nreads)];


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

testLikeRealData_genericChooseCorrection(outputName,['toyExample/sim1/reads/',basicDirNameForFile],auxDataTest)

unix(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_',outputName,'/run_sol_',outputName,'_noCorrection_',num2str(auxDataTest.readLength)])

unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/12Samples/Solexa/results/',auxDataTest.basicDirNameForFile])
unix(['cp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/sol_',outputName{i},'/sol_',outputName{i},'_noCorrection_',num2str(auxDataTest.readLength),'/sol_',outputName{i},'_noCorrection_',num2str(auxDataTest.readLength),'.mat',' /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/12Samples/Solexa/results/',auxDataTest.basicDirNameForFile])

unix(['cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/12Samples/Solexa/results/',';tar cf results_', ...
      auxDataTest.basicDirNameForFile,'.tar ',auxDataTest.basicDirNameForFile])


%%%%%%%%
% run 350 no noise
clear
if  ~isempty(findstr(getenv('HOSTNAME'),'openu'))
  userDir = getuserdir;
else
  userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
end

readLength = 350;
addNoiseFlag = 0;
Nreads = 10^4;
dataName = ['reads350_readlen_',num2str(readLength),'_noise_',num2str(addNoiseFlag),'_Nreads_',num2str(Nreads)];

auxData = struct;
auxData.dataName = dataName;
auxData.userDir = userDir;
auxData.basicDir = [auxData.userDir,'/'];
auxData.basicName = ['test_',dataName];
auxData.databaseName = [auxData.userDir,'/CS/BAC/12Samples/454/idx/old_fasta_files_without_ambiguous/454idx.fa'];
auxData.blastPath = [auxData.userDir,'/CS/BAC/12Samples/454/idx/blast-2.2.17/bin/'];
auxData.dataDir = [auxData.userDir,'/CS/BAC/toyExample/sim1/reads'];
auxData.numProcessors = 2;
auxData.queueName = 'hour';

unix(['mkdir ',auxData.userDir,'/CS/BAC/alignToDatabase/tmpRuns/alignTmp_',auxData.dataName])
unix(['mkdir ',auxData.userDir,'/CS/BAC/alignToDatabase/dataForSim/alignTmp_',auxData.dataName])
unix(['mkdir ',auxData.basicDir,';mkdir ',auxData.basicDir,'/tmp'])

alignToDatabse(auxData)

% end run 350


% end run no noise case

%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%

% run with noise


%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%

% run 350 with noise
clear
if  ~isempty(findstr(getenv('HOSTNAME'),'openu'))
  userDir = getuserdir;
else
  userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
end

readLength = 350;
addNoiseFlag = 1;
Nreads = 10^4;
dataName = ['reads350_readlen_',num2str(readLength),'_noise_',num2str(addNoiseFlag),'_Nreads_',num2str(Nreads)];

auxData = struct;
auxData.dataName = dataName;
auxData.userDir = userDir;
auxData.basicDir = [auxData.userDir,'/'];
auxData.basicName = ['test_',dataName];
auxData.databaseName = [auxData.userDir,'/CS/BAC/12Samples/454/idx/old_fasta_files_without_ambiguous/454idx.fa'];
auxData.blastPath = [auxData.userDir,'/CS/BAC/12Samples/454/idx/blast-2.2.17/bin/'];
auxData.dataDir = [auxData.userDir,'/CS/BAC/toyExample/sim1/reads'];
auxData.numProcessors = 100;
auxData.queueName = 'hour';

unix(['mkdir ',auxData.userDir,'/CS/BAC/alignToDatabase/tmpRuns/alignTmp_',auxData.dataName])
unix(['mkdir ',auxData.userDir,'/CS/BAC/alignToDatabase/dataForSim/alignTmp_',auxData.dataName])
%unix(['mkdir ',auxData.basicDir,';mkdir ',auxData.basicDir,'/tmp'])

alignToDatabse(auxData)




% run 750 with noise case - no correction
clear
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

dataName = ['reads750_readlen_',num2str(readLength),'_noise_',num2str(addNoiseFlag),'_Nreads_',num2str(Nreads)];

basicDirNameForFile = [dataName,'_correction_',num2str(auxDataTest.correctMeasurementForReadErrorsFlag)]
outputName = ['test_',basicDirNameForFile];

testLikeRealData_genericChooseCorrection(outputName,['toyExample/sim1/reads/',dataName],auxDataTest)

unix(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_',outputName,'/run_sol_',outputName,'_noCorrection_readlen_',num2str(auxDataTest.readLength)])

% run 750 with noise case - with correction
clear
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

dataName = ['reads750_readlen_',num2str(readLength),'_noise_',num2str(addNoiseFlag),'_Nreads_',num2str(Nreads)];

basicDirNameForFile = [dataName,'_correction_',num2str(auxDataTest.correctMeasurementForReadErrorsFlag)]
outputName = ['test_',basicDirNameForFile];

testLikeRealData_genericChooseCorrection(outputName,['toyExample/sim1/reads/',dataName],auxDataTest)

unix(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_',outputName,'/run_sol_',outputName,'_withCorrection_readlen_',num2str(auxDataTest.readLength)])

% copy the results to one directory
unix(['cp '])


% more:
% run with noise
% high number of reads, low number of reads
% lower number of reads in 454


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create tree
clear 

if  ~isempty(findstr(getenv('HOSTNAME'),'openu'))
  userDir = getuserdir;
else
  userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
end


load ~/CS/BAC/toyExample/sim1/forTree/header_sim1_toy

Header750 = load('~/CS/BAC/primers750_primer_tail/bac16s_primers750_primer_tail_full_without_ambiguous');

clear posIn750 currseq lenSeq
for i=1:length(header_sim_toy)
  i
  for j=1:length(Header750.Header_750_tail)
    for k=1:length(Header750.Header_750_tail{j})
      if ~isempty(findstr(Header750.Header_750_tail{j}{k},header_sim_toy{i}{1}))
        posIn750{i} = j;
        currseq{i} = Header750.Sequence_750_tail{j};
        lenSeq(i) = length(currseq{i});
      end
    end
  end
end



reads_uni = char(55*ones(length(currseq),max(lenSeq)));
for i=1:length(lenSeq)
  reads_uni(i,1:lenSeq(i)) = currseq{i};
end

save ~/CS/BAC/toyExample/sim1/forTree/seqs_header_sim1_toy reads_uni lenSeq

auxData = struct;
auxData.userDir = userDir;
auxData.blastPath = [auxData.userDir,'/CS/BAC/12Samples/454/idx/blast-2.2.17/bin/'];
dataName = 'seqs_header_sim1_toy';
auxData.basicName = ['test_',dataName];
auxData.basicDir = [auxData.userDir,'/'];
auxData.databaseName = [auxData.userDir,'/CS/BAC/primers750_primer_tail/fastaFile/bac16s_primers750_primer_tail_full_without_ambiguous.fa'];
auxData.blastPath = [auxData.userDir,'/CS/BAC/12Samples/454/idx/blast-2.2.17/bin/'];
auxData.saveFileName = dataName;
auxData.dataDir = [auxData.userDir,'/CS/BAC/toyExample/sim1/forTree'];
auxData.dataName = dataName;
auxData.numProcessors = 200;
auxData.queueName = 'hour';

%unix(['mkdir ',auxData.userDir,'/CS/BAC/alignToDatabase/tmpRuns/alignTmp_',auxData.dataName])
%unix(['mkdir ',auxData.userDir,'/CS/BAC/alignToDatabase/dataForSim/alignTmp_',auxData.dataName])
%unix(['mkdir ',auxData.basicDir,';mkdir ',auxData.basicDir,'/tmp'])

%alignToDatabseWithThreshold(auxData)

myFindSeqInDB_genericWithThreshold(auxData.basicName,auxData.basicDir,auxData.databaseName,auxData.blastPath,1,length(currseq),auxData.saveFileName,auxData.dataDir,auxData.dataName)


% save the sequences and Headers
res = load(auxData.saveFileName);
clear HeadAndSeq_res
for i=1:length(res.idx)
  for j=1:length(res.idx{i})
    HeadAndSeq_res{i}.Header{j} = Header750.Header_750_tail{res.idx{i}(j)};
    HeadAndSeq_res{i}.Sequence{j} = Header750.Sequence_750_tail{res.idx{i}(j)};
    HeadAndSeq_res{i}.idx = res.idx{i};
    HeadAndSeq_res{i}.dist = res.dist{i};
  end
end

save([auxData.dataDir,'/HeaderAndSeq_sim1'],'HeadAndSeq_res')

