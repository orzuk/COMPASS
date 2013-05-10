clear

[Header_R1, Sequence_R1] = fastqread('~/CS/BAC/CMRNA5/data/CMRNA5_CTTGTA_R1.fastq');
[Header_R2, Sequence_R2] = fastqread('~/CS/BAC/CMRNA5/data/CMRNA5_CTTGTA_R2.fastq');


[reads_R1,non] = extractReadsPeterFiles(Sequence_R1);
reads_R1(non) = [];

[reads_R2,non] = extractReadsPeterFiles(Sequence_R2);
reads_R2(non) = [];

reads =[reads_R1,reads_R2];







%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find uniqueReads
m4sort=sort(reads);
[uni_reads,ia]=unique(m4sort,'first');
[~,ib]=unique(m4sort,'last');
uniqueReads_length = ib-ia+1;

uniqueReads = pack_seqs(uni_reads,64);
uniqueReads = cell2mat(uniqueReads);

tmp_uni_reads = char(zeros(length(uni_reads),length(uni_reads{1})));
for i=1:length(uni_reads)
  tmp_uni_reads(i,:) = uni_reads{i};
end
uni_reads = tmp_uni_reads;

save ~/CS/BAC/CMRNA5/data/dataCMRNA5 uni_reads uniqueReads_length uniqueReads

%uni_reads = uni_reads(1:10000,:);
%uniqueReads = uniqueReads(1:10000,:);
%uniqueReads_length = uniqueReads_length(1:10000);
%save ~/CS/BAC/CMRNA5/data/dataCMRNA5_part uni_reads uniqueReads uniqueReads_length



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% align the reads
clear
%userDir = getuserdir;
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';

sampleName = 'CMRNA5';

auxData = struct;
auxData.queueName = 'hour';
auxData.readLength = 101;
auxData.searchForwardAndReverseFlag = 0; %genric - pay attanetion if read are from both starnads!!!  
auxData.readsFileName = ['dataCMRNA5'];
auxData.inputDirName = [userDir,'/CS/BAC/CMRNA5/data'];
auxData.tmpBasicFileName = ['alignTmp_',num2str(auxData.readLength),'_',sampleName];


ul = load([auxData.inputDirName,'/',auxData.readsFileName],'uniqueReads_length');
n = length(ul.uniqueReads_length);

auxData.n = n;
auxData.numProcessors = round(n/1000);
auxData.numProcessors

auxData.charVectorFile = [userDir,'/CS/BAC/12Samples/Solexa/data/forAlign_full16S/dataForAlign_full16S'];


unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/findAlignment_generic/tmpRuns/alignTmp_',num2str(auxData.readLength),'_',sampleName])
unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/findAlignment_generic/dataForSim/alignTmp_',num2str(auxData.readLength),'_',sampleName])



findPositionOfReads_generic(auxData);
  

% filter those that are not 16S
clear
load ~/CS/BAC/CMRNA5/data/dataCMRNA5
load('~/CS/BAC/CMRNA5/data/alignTmp_101_CMRNA5_finalRes')

a = find(POS);
length(a)


b = unique(POS(a));
num = zeros(length(b),1);
for i=1:length(b)
  c = find(POS(a)==b(i));
  num(i) = sum(uniqueReads_length(a(c)));
end


z = round(b);
u = unique(z);
for i=1:length(u)
  n(i) = sum(num(find(z==u(i))));
end

plot(u,n,'.','markersize',10)
set(gca,'fontweight','bold','fontsize',15)
title('coverage along the gene')
print -dpdf ~/CS/BAC/CMRNA5/coverage_along_gene

% plot the distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare the reads - only aligned - run without correction
if 1==2
  uniqueReads = uniqueReads(a,:);
  uniqueReads_length = uniqueReads_length(a);
  uni_reads = uni_reads(a,:);
  save ~/CS/BAC/CMRNA5/data/uniqueAlignedReadsCMRNA5 uniqueReads uniqueReads_length uni_reads
  clear uniqueReads uniqueReads_length uniqueReads uni_reads
end

% prepare the 90 version


% run unnormalized
runDataCMRNA5_no_eq_no_correction('testCMRNA5_alignedReads_noNor')

% run normalized

% run 90 without normalization

% run 90 with normalization



%two runs

clear
load ~/CS/BAC/sol_testCMRNA5_alignedReads_testCMRNA5_alignedReads_noCorrection_101
load ~/CS/BAC/full16S/bac16s_full_without_ambiguous
load ~/CS/BAC/CMRNA5/uniqueAlignedReadsCMRNA5

[freq,ind] = sort(abs(found{end}),'descend');

for i=1:length(freq)
  fprintf('%1.2f',freq(i));
  fprintf('%s\n',Header_uni{ind(i)});
  pause
end


a = find(abs(found{end})>10^-3);


% test l1
userDir = getuserdir;
basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
basicSeqKey = [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
readLength = 101;
tmpInd = find(abs(found{end})>10^-3);
[normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(readLength,tmpInd,basicSeqNameDir,basicSeqKey);
  

length(find(b==1))

b = sum(normalizedBac,2);
b = find(b==1);
[x,y] = find(normalizedBac(b,:)==1);


[vals inds num_dups] = get_duplicates2(y);

plot(num_dups,'.')
title('number of unique kmers of each bateria')


s = zeros(1,301);
for i=1:length(vals)
  s(vals(i)) = sum(fracRelevantReads(b(inds{i})));
end

plot(s,'.')
title('number of actual unique reads for each bacteria')


dataIn = struct;
[fracRelevantReads,sumRelevantReads] = currReadsFourth(uniqueReads,uniqueReads_length,values,0,dataIn);

[x1_2_nor] = testL1_2(normalizedBac./max(fracRelevantReads),fracRelevantReads./max(fracRelevantReads));    
[x1_2_no_nor] = testL1_2(normalizedBac,fracRelevantReads);    


figure(2)
subplot(1,3,1)
plot(abs(found{end}(a)),'.')
subplot(1,3,2)
plot(x1_2_nor,'.');title('l1 nor')
subplot(1,3,3)
plot(x1_2_no_nor,'.');title('l1 no nor')


save ~/CS/BAC/CMRNA5/dataForComparison tmpInd found x1_2_nor  normalizedBac fracRelevantReads values

% freq_l2 = abs(found{end}(tmpInd));save ~/CS/BAC/CMRNA5/dataL1_L2_CMRNA5 freq_l2 normalizedBac fracRelevantReads x1_2_nor x1_2_no_nor

% distance matrix
load(basicSeqKey)
[HeaderAll,SequenceAll] = loadSeq(tmpInd,basicSeqNameDir,basicSeqKey);
ham = zeros*ones(length(tmpInd));
for i=1:length(tmpInd)-1
  i
  seq_tmpInd = int2nt(unpack_seqs(SequenceAll{i},len_uni(tmpInd(i)),64));
  for j=i+1:length(tmpInd)
    seq_rec = int2nt(unpack_seqs(SequenceAll{j},len_uni(tmpInd(j)),64));
    [cscore,algn]=swalign(seq_tmpInd,seq_rec,'Alphabet','NT');
    ham(i,j) = length(find(algn=='|'));
    
    %ham(i,j)
    %pause
  end
  length(find((ham)))
end

ham = ham+ham';

imagesc(ham)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
% simulated data
clear
fileName = 'data_reads_peter_29';
inputName = 'reads-peter-29';
outputName = 'test_peter_29';

readLength = 100;
createUniqueReadsSimFollowingCMRNA5(['~/CS/BAC/CMRNA5/',inputName],['~/CS/BAC/CMRNA5/',fileName],readLength)

% copy to br
testSimFromReadsLikeRealData_fullDatabse(outputName,['CMRNA5/',fileName],50,100)

unix(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_',outputName,'_',outputName,'/run_sol_',outputName,'_',outputName,'_noCorrection_100'])

unix(['cp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/sol_',outputName,'_',outputName,'/sol_',outputName,'_',outputName,'_noCorrection_100/sol_',outputName,'_',outputName,'_noCorrection_100.mat',' /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/CMRNA5/'])

extractResSimCMRNA5(outputName,inputName,fileName,100)

%%%%%%%%%%%%%%
clear
fileName = 'data_reads_peter_29_2';
inputName = 'reads-peter-29-2';
outputName = 'test_peter_29-2';

readLength = 100;
createUniqueReadsSimFollowingCMRNA5(['~/CS/BAC/CMRNA5/',inputName],['~/CS/BAC/CMRNA5/',fileName],readLength)

% copy to br
testSimFromReadsLikeRealData_fullDatabse(outputName,['CMRNA5/',fileName],50,100)

unix(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_',outputName,'_',outputName,'/run_sol_',outputName,'_',outputName,'_noCorrection_100'])

unix(['cp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/sol_',outputName,'_',outputName,'/sol_',outputName,'_',outputName,'_noCorrection_100/sol_',outputName,'_',outputName,'_noCorrection_100.mat',' /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/CMRNA5/'])

extractResSimCMRNA5(outputName,inputName,fileName,100)


%%%%%%%%%%%%%%%%%
clear
fileName = 'data_reads_peter_301';
inputName = 'reads-peter-301';
outputName = 'test_peter_301';

readLength = 100;
createUniqueReadsSimFollowingCMRNA5(['~/CS/BAC/CMRNA5/',inputName],['~/CS/BAC/CMRNA5/',fileName],readLength)

% copy to br
testSimFromReadsLikeRealData_fullDatabse(outputName,['CMRNA5/',fileName],50,100)

unix(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_',outputName,'_',outputName,'/run_sol_',outputName,'_',outputName,'_noCorrection_100'])

unix(['cp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/sol_',outputName,'_',outputName,'/sol_',outputName,'_',outputName,'_noCorrection_100/sol_',outputName,'_',outputName,'_noCorrection_100.mat',' /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/CMRNA5/'])

extractResSimCMRNA5(outputName,inputName,fileName,100)

%%%%%%%%%%%%%%%
clear
fileName = 'data_reads_peter_301_2';
inputName = 'reads-peter-301-2';
outputName = 'test_peter_301-2';

readLength = 100;
createUniqueReadsSimFollowingCMRNA5(['~/CS/BAC/CMRNA5/',inputName],['~/CS/BAC/CMRNA5/',fileName],% copy to br
testSimFromReadsLikeRealData_fullDatabse(outputName,['CMRNA5/',fileName],50,100)

unix(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_',outputName,'_',outputName,'/run_sol_',outputName,'_',outputName,'_noCorrection_100'])

unix(['cp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/sol_',outputName,'_',outputName,'/sol_',outputName,'_',outputName,'_noCorrection_100/sol_',outputName,'_',outputName,'_noCorrection_100.mat',' /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/CMRNA5/'])

extractResSimCMRNA5(outputName,inputName,fileName,100)


%%%%%%%%%%%%%%%%%%%%%%%%%%%55

clear
fileName = {'data_reads_peter_301_noise_0','data_reads_peter_301_noise_1','data_reads_peter_301_noise_100','data_reads_peter_29_noise_0','data_reads_peter_29_noise_1','data_reads_peter_29_noise_100'};
inputName = {'reads-peter-301-noise-0','reads-peter-301-noise-1','reads-peter-301-noise-100','reads-peter-29-noise-0','reads-peter-29-noise-1','reads-peter-29-noise-100'};
outputName = {'test_peter_301-noise-0','test_peter_301-noise-1','test_peter_301-noise-100','test_peter_29-noise-0','test_peter_29-noise-1','test_peter_29-noise-100'};



%fileName = {'data_reads_peter_29_noise_0','data_reads_peter_29_noise_1','data_reads_peter_29_noise_100'};
%inputName = {'reads-peter-29-noise-0','reads-peter-29-noise-1','reads-peter-29-noise-100'};
%outputName = {'test_peter_29-noise-0','test_peter_29-noise-1','test_peter_29-noise-100'};


readLength = 100;

for i=1:length(fileName)
  createUniqueReadsSimFollowingCMRNA5(['~/CS/BAC/CMRNA5/',inputName{i}],['~/CS/BAC/CMRNA5/',fileName{i}],readLength);
end

% move to br
for i=1:length(fileName)
  testSimFromReadsLikeRealData_fullDatabse(outputName{i},['CMRNA5/',fileName{i}],50,100)
  unix(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_',outputName{i},'_',outputName{i},'/run_sol_',outputName{i},'_',outputName{i},'_noCorrection_100'])
  pause(10)
end

for i=1:length(fileName)
  unix(['cp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/sol_',outputName{i},'_',outputName{i},'/sol_',outputName{i},'_',outputName{i},'_noCorrection_100/sol_',outputName{i},'_',outputName{i},'_noCorrection_100.mat',' /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/CMRNA5/'])

end

for i=1:length(fileName)
  extractResSimCMRNA5(outputName{i},inputName{i},fileName{i},100)
  %pause
end




%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
% send to CMRNA5

clear
load ~/CS/BAC/sol_testCMRNA5_alignedReads_testCMRNA5_alignedReads_noCorrection_101
load ~/CS/BAC/full16S/bac16s_full_without_ambiguous
load ~/CS/BAC/CMRNA5/uniqueAlignedReadsCMRNA5


userDir = getuserdir;
basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];
basicSeqKey = [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
readLength = 101;
tmpInd = find(abs(found{end})>10^-3);
[normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(readLength,tmpInd,basicSeqNameDir,basicSeqKey);
  

dataIn = struct;
[fracRelevantReads,sumRelevantReads] = currReadsFourth(uniqueReads,uniqueReads_length,values,0,dataIn);

[x1_2_nor] = testL1_2(normalizedBac./max(fracRelevantReads),fracRelevantReads./max(fracRelevantReads));    

a = find(x1_2_nor>10^-2);
[junk,ind] = sort(x1_2_nor(a),'descend');
a = a(ind);

curr_ind = tmpInd(a);


fid = fopen('~/CS/BAC/CMRNA5/results_4_4_12.xls','w');
fprintf(fid,'freq\t Header\n');
for i=1:length(ind)
  fprintf(fid,'%1.2f\t %s\n',junk(i),Header_uni{curr_ind(i)});
end
fclose(fid)


Header_CMRNA5 = Header_uni(curr_ind);
Sequence_CMRNA5 = Sequence_uni(curr_ind);
freq = junk;

save ~/CS/BAC/CMRNA5/res_4_4_12 Header_CMRNA5 Sequence_CMRNA5 freq

fastawrite('~/CS/BAC/CMRNA5/res_4_4_12.fa',Header_CMRNA5,Sequence_CMRNA5);

if 1==2
  clear
  load ~/CS/BAC/sol_testCMRNA5_alignedReads_testCMRNA5_alignedReads_noCorrection_101
  load ~/CS/BAC/full16S/bac16s_full_without_ambiguous
  tmpInd = find(abs(found{end})>10^-3);
  
  freq301 = abs(found{end}(tmpInd));
  
  Header_301 = Header_uni(tmpInd);
  Sequence_301 = Sequence_uni(tmpInd);
  fastawrite('~/CS/BAC/CMRNA5/res_CMRNA5_301.fa',Header_301,Sequence_301);

  save ~/CS/BAC/CMRNA5/res_301 Header_301 Sequence_301 freq301
  
end


used 101*2 reads, 1000000, 300000 unique- 


s = sum(normalizedBac,2);


s1 = find(s==1);
s2 = find(s==2);

for i=1:size(normalizedBac,2)
  a = find(normalizedBac(:,i));
  b = intersect(a,s1);
  sum_zero(i) = length(find(fracRelevantReads(b)==0));
  sum_non_zero(i) = sum(fracRelevantReads(b));

  sum_total(i) = sum(fracRelevantReads(a));
  
end


plot(sum_total,sum_zero,'.')
xlabel('total reads')
ylabel('num kmers with zero reads')

figure
plot(sum_total,sum_non_zero,'.')
xlabel('total reads')
ylabel('unique reads')

figure
plot(sum_total,sum_non_zero,'.')
hold on
c = find(sum_non_zero>50);
plot(sum_total(c),sum_non_zero(c),'ro')















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%5555
% create fasta files of reads themselves
clear
load ~/CS/BAC/CMRNA5/uniqueAlignedReadsCMRNA5


reads = int2nt(unpack_seqs(uniqueReads,101*ones(size(uniqueReads,1),1),64));

r = cell(size(reads,1),1);
h = r;
for i=1:size(reads,1)
  r{i} = reads(i,:);
  h{i} = num2str(i);
end

delete('~/CS/BAC/CMRNA5/sampleReads/sample1_1000.fa')

fastawrite('~/CS/BAC/CMRNA5/sampleReads/sample1_10000.fa',h(1:10000),r(1:10000));
fastawrite('~/CS/BAC/CMRNA5/sampleReads/sample10001_20000.fa',h(10001:20000),r(10001:20000));
fastawrite('~/CS/BAC/CMRNA5/sampleReads/sample20001_31575.fa',h(20001:31575),r(20001:31575));


%%%%%%%%%%%%%%%%
clear val fam

num_unique_reads = size(uniqueReads,1);
fileName = ['~/CS/BAC/CMRNA5/sampleReads/sample1_10000.fa_download.txt']
[fam_tmp,val_tmp]=readClassifedReads(fileName,num_unique_reads);
fam(1:10000) = fam_tmp(1:10000);
val(1:10000) = val_tmp(1:10000);

fileName = ['~/CS/BAC/CMRNA5/sampleReads/sample10001_20000.fa_download.txt']
[fam_tmp,val_tmp]=readClassifedReads(fileName,num_unique_reads);
fam(10001:20000) = fam_tmp(10001:20000);
val(10001:20000) = val_tmp(10001:20000);

fileName = ['~/CS/BAC/CMRNA5/sampleReads/sample20001_31575.fa_download.txt']
[fam_tmp,val_tmp]=readClassifedReads(fileName,num_unique_reads);
fam(20001:31575) = fam_tmp(20001:31575);
val(20001:31575) = val_tmp(20001:31575);


u = unique(fam);

for i=1:length(u)
  a = find(ismember(fam, u{i})==1);
  sum_reads(i) = sum(uniqueReads_length(a));
  numRepeats(i) = sum(uniqueReads_length(a)'.*val(a)/100);
end

save ~/CS/BAC/CMRNA5/sampleReads/rawDataClassification sum_reads numRepeats u fam val



plot(numRepeats,'.')
set(gca,'fontweight','bold','fontsize',15,'xtick',1:length(u),'xticklabel',u','position',[0.1 0.25 0.8 0.65])
rotateticklabel(gca,45);
ylabel('#reads')
print -dpdf ~/CS/BAC/CMRNA5/sampleReads/RDPbasedClassification


plot(sum_reads,'.')
set(gca,'fontweight','bold','fontsize',15,'xtick',1:length(u),'xticklabel',u','position',[0.1 0.25 0.8 0.65])
rotateticklabel(gca,45);
ylabel('#reads')
print -dpdf ~/CS/BAC/CMRNA5/sampleReads/RDPbasedClassification_sum_appearances


res_21 = load('~/CS/BAC/CMRNA5/res_4_4_12'); 
fileName = ['~/CS/BAC/CMRNA5/res_4_4_12.fa_download.txt']
[fam_21,val_21]=readClassifedReads(fileName,21,'roling');


res_301 = load('~/CS/BAC/CMRNA5/res_301');
fileName = ['~/CS/BAC/CMRNA5/res_CMRNA5_301.fa_download.txt'];
[fam_301,val_301]=readClassifedReads(fileName,301,'roling');


for i=1:length(u)
  a = find(ismember(fam_21, u{i})==1);
  freq_21(i) = sum(res_21.freq(a));
  
  a = find(ismember(fam_301, u{i})==1);
  freq_301(i) = sum(res_301.freq301(a));
  
end

figure(4);clf
subplot(3,1,1)
plot(freq_21,'o','markersize',6,'markerfacecolor','b')
set(gca,'fontweight','bold','fontsize',15,'ytick',[0:0.1:1],'ylim',[0 1])
ylabel('phyla freq')
title('version 1')


subplot(3,1,2)
plot(freq_301,'o','markersize',6,'markerfacecolor','b')
set(gca,'fontweight','bold','fontsize',15,'ytick',[0:0.1:1],'ylim',[0 1])
ylabel('phyla freq')
title('version 2')


subplot(3,1,3)
plot(sum_reads,'o','markersize',6,'markerfacecolor','b')
set(gca,'fontweight','bold','fontsize',15,'xtick',1:length(u),'xticklabel',u')
ylabel('#reads')
title('classification of raw reads')
%set(gca,'fontweight','bold','fontsize',15,'xtick',1:length(u),'xticklabel',u','position',[0.1300    0.7093    0.7750    0.2068])
z = rotateticklabel(gca,45);
set(z,'fontsize',10)


set(gcf,'paperposition',[0 1 9, 10])
print -dpdf ~/CS/BAC/CMRNA5/sampleReads/phyla_classification_comparison


%%%%%%%%%%%%%%%%%%%%%%5555
%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%5555
% 24.4

%%%%%%%%%%%%%%%%%
% test normalized reads in CMRNA5 
clear
fileName = 'myFormat_normalizedreads';
inputName = 'normalizedreads';
outputName = 'test_normalizedreads';
readLength = 101;

disp('based on all reads - with those that are not aligned, but does not correct for that')


origData = load('~/CS/BAC/CMRNA5/dataCMRNA5') 
load ~/CS/BAC/CMRNA5/normalizedreads

a = find(normreads);

uniqueReads_length = normreads(a);
uniqueReads = origData.uniqueReads(a,:);

save(['~/CS/BAC/CMRNA5/',fileName],'uniqueReads','uniqueReads_length')


% copy to br
testSimFromReadsLikeRealData_fullDatabse(outputName,['CMRNA5/',fileName],50,readLength)

unix(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_',outputName,'_',outputName,'/run_sol_',outputName,'_',outputName,'_noCorrection_',num2str(readLength)])

unix(['cp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/sol_',outputName,'_',outputName,'/sol_',outputName,'_',outputName,'_noCorrection_',num2str(readLength),'/sol_',outputName,'_',outputName,'_noCorrection_',num2str(readLength),'.mat',' /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/CMRNA5/'])

extractResSimCMRNA5(outputName,inputName,fileName,readLength)

%%%%%%%%%%%%%%%%%%%
clear
load ~/CS/BAC/CMRNA5/simRes/res_normalizedreads
load ~/CS/BAC/full16S/bac16s_full_without_ambiguous


a = find(x1_2_nor>10^-3);
[junk,ind] = sort(x1_2_nor(a),'descend');
a = a(ind);

curr_ind = tmpInd(a);


Header_CMRNA5 = Header_uni(curr_ind);
Sequence_CMRNA5 = Sequence_uni(curr_ind);
freq = junk;



fileName = '~/CS/BAC/CMRNA5/res_25_4_12.fa';
delete(fileName)
fastawrite(fileName,Header_CMRNA5,Sequence_CMRNA5);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do CMRNA5 with 75

clear
alignReads = load('~/CS/BAC/CMRNA5/uniqueAlignedReadsCMRNA5'); 
readLengthPost = 75;
readLengthPre = 101;
fileName = 'myFormat_CMRNA575';
inputName = 'CMRNA575';
outputName = 'test_CMRNA575';


r = int2nt(unpack_seqs(alignReads.uniqueReads,101,64));
uni_reads = cell(size(r,1),1);
for i=1:size(r,1)
  uni_reads{i} = r(i,:);
end



r_post = repmat(char(zeros(size(uni_reads,1),size(r,2))),readLengthPre-readLengthPost+1,1);
r_post = r_post(:,1:readLengthPost);
length_r_post = zeros(size(r_post,1),1);
k = 1;
for i=1:size(uni_reads,1)
  if mod(i,10^5)==1
    i
  end
  for j=1:readLengthPre-readLengthPost+1
    r_post(k,:) = uni_reads{i}(j:j+readLengthPost-1);
    length_r_post(k) = alignReads.uniqueReads_length(i);
    k = k+1;
  end
end

r_post(k+1:end,:) = [];


readChar = cell(size(r_post,1),1);
for i=1:size(r_post,1)
  readChar{i} = r_post(i,:);
end


[vals1 inds1 num_dups] = get_duplicates2(readChar);

uniqueReads = pack_seqs(vals1,64);
uniqueReads = cell2mat(uniqueReads);

uniqueReads_length = zeros(size(vals1,1),1);
for i=1:size(uniqueReads_length,1)
  uniqueReads_length(i) = sum(length_r_post(inds1{i}));
end


save(['~/CS/BAC/CMRNA5/',fileName],'uniqueReads','uniqueReads_length')



% copy to br
testSimFromReadsLikeRealData_fullDatabse(outputName,['CMRNA5/',fileName],50,readLengthPost)

unix(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_',outputName,'_',outputName,'/run_sol_',outputName,'_',outputName,'_noCorrection_',num2str(readLengthPost)])

unix(['cp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/sol_',outputName,'_',outputName,'/sol_',outputName,'_',outputName,'_noCorrection_',num2str(readLengthPost),'/sol_',outputName,'_',outputName,'_noCorrection_',num2str(readLengthPost),'.mat',' /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/CMRNA5/'])

extractResSimCMRNA5(outputName,inputName,fileName,readLengthPost)

%%%%%%%%%%%%%%
clear
load ~/CS/BAC/CMRNA5/simRes/res_CMRNA575
load ~/CS/BAC/full16S/bac16s_full_without_ambiguous


a = find(x1_2_nor>10^-3);
[junk,ind] = sort(x1_2_nor(a),'descend');
a = a(ind);

curr_ind = tmpInd(a);


Header_CMRNA5 = Header_uni(curr_ind);
Sequence_CMRNA5 = Sequence_uni(curr_ind);
freq = junk;



fileName = '~/CS/BAC/CMRNA5/res_CMRNA575.fa';
delete(fileName)
fastawrite(fileName,Header_CMRNA5,Sequence_CMRNA5);

%%%%%%%%%%%%%%%%%
% do with 90
clear
alignReads = load('~/CS/BAC/CMRNA5/uniqueAlignedReadsCMRNA5'); 
readLengthPost = 90;
readLengthPre = 101;
fileName = 'myFormat_CMRNA590';
inputName = 'CMRNA590';
outputName = 'test_CMRNA590';


r = int2nt(unpack_seqs(alignReads.uniqueReads,101,64));
uni_reads = cell(size(r,1),1);
for i=1:size(r,1)
  uni_reads{i} = r(i,:);
end



r_post = repmat(char(zeros(size(uni_reads,1),size(r,2))),readLengthPre-readLengthPost+1,1);
r_post = r_post(:,1:readLengthPost);
length_r_post = zeros(size(r_post,1),1);
k = 1;
for i=1:size(uni_reads,1)
  if mod(i,10^5)==1
    i
  end
  for j=1:readLengthPre-readLengthPost+1
    r_post(k,:) = uni_reads{i}(j:j+readLengthPost-1);
    length_r_post(k) = alignReads.uniqueReads_length(i);
    k = k+1;
  end
end

r_post(k+1:end,:) = [];


readChar = cell(size(r_post,1),1);
for i=1:size(r_post,1)
  readChar{i} = r_post(i,:);
end


[vals1 inds1 num_dups] = get_duplicates2(readChar);

uniqueReads = pack_seqs(vals1,64);
uniqueReads = cell2mat(uniqueReads);

uniqueReads_length = zeros(size(vals1,1),1);
for i=1:size(uniqueReads_length,1)
  uniqueReads_length(i) = sum(length_r_post(inds1{i}));
end


save(['~/CS/BAC/CMRNA5/',fileName],'uniqueReads','uniqueReads_length')


% copy to br
testSimFromReadsLikeRealData_fullDatabse(outputName,['CMRNA5/',fileName],50,readLengthPost)

unix(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/sol_',outputName,'_',outputName,'/run_sol_',outputName,'_',outputName,'_noCorrection_',num2str(readLengthPost)])

unix(['cp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/sol_',outputName,'_',outputName,'/sol_',outputName,'_',outputName,'_noCorrection_',num2str(readLengthPost),'/sol_',outputName,'_',outputName,'_noCorrection_',num2str(readLengthPost),'.mat',' /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/CMRNA5/'])

extractResSimCMRNA5(outputName,inputName,fileName,readLengthPost)








%%%%%%%%%%%%%%%55
% 25.4.12
% run with correction
runDataCMRNA5_no_eq_with_correction('testCMRNA5_withCorrection','dataCMRNA5')

% did not work

%%%%%%%%%%%%%%%%%

% prepare reads from both legs and look at distribution
clear

[Header, Sequence] = fastaread('~/CS/BAC/CMRNA5/A_Amp_16S.fna');

L = zeros(length(Sequence),1);
for i=1:length(Sequence)
  L(i) = length(Sequence{i});
end

I = find(L==202);
reads_right = cell(length(I),1);
reads_left = reads_right;
for i=1:length(I)
  reads_right{i,:} = Sequence{I(i)}(1:101);
  reads_left{i,:} =  Sequence{I(i)}(102:202);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find uniqueReads
m4sort=sort(reads_right);
[uni_reads,ia]=unique(m4sort,'first');
[~,ib]=unique(m4sort,'last');
uniqueReads_length = ib-ia+1;

uniqueReads = pack_seqs(uni_reads,64);
uniqueReads = cell2mat(uniqueReads);

save ~/CS/BAC/CMRNA5/dataCMRNA5_right uni_reads uniqueReads_length uniqueReads

% find uniqueReads
m4sort=sort(reads_left);
[uni_reads,ia]=unique(m4sort,'first');
[~,ib]=unique(m4sort,'last');
uniqueReads_length = ib-ia+1;

uniqueReads = pack_seqs(uni_reads,64);
uniqueReads = cell2mat(uniqueReads);

save ~/CS/BAC/CMRNA5/dataCMRNA5_left uni_reads uniqueReads_length uniqueReads


%%%%%%%%
clear 
right = load('~/CS/BAC/CMRNA5/dataCMRNA5_right');
left = load('~/CS/BAC/CMRNA5/dataCMRNA5_left');


%%%%%%%%
% in br

% align the reads - right
clear
auxData = struct;
auxData.queueName = 'hour';
auxData.readLength = 101;

n = 115350;
auxData.n = n;
auxData.numProcessors = round(n/1000);
auxData.numProcessors


auxData.dataCMRNA5Name = 'dataCMRNA5_right';

sampleName = 'CMRNA5_right';
unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/findAlignment_CMRNA5/tmpRuns/alignTmp_',num2str(auxData.readLength),'_',sampleName])
unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/findAlignment_CMRNA5/dataForSim/alignTmp_',num2str(auxData.readLength),'_',sampleName])

findPositionOfReads_CMRNA5(sampleName,auxData);

% align the reads - left
clear
auxData = struct; 
auxData.queueName = 'hour';
auxData.readLength = 101;

n = 121384;
auxData.n = n;
auxData.numProcessors = round(n/1000);
auxData.numProcessors

auxData.dataCMRNA5Name = 'dataCMRNA5_left';

sampleName = 'CMRNA5_left';
unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/findAlignment_CMRNA5/tmpRuns/alignTmp_',num2str(auxData.readLength),'_',sampleName])
unix(['mkdir /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/findAlignment_CMRNA5/dataForSim/alignTmp_',num2str(auxData.readLength),'_',sampleName])

findPositionOfReads_CMRNA5(sampleName,auxData);

%%%%%%%%% right
clear
right = load('~/CS/BAC/CMRNA5/dataCMRNA5_right') 
right_res = load('~/CS/BAC/CMRNA5/alignTmp_101_CMRNA5_right_finalRes');

left = load('~/CS/BAC/CMRNA5/dataCMRNA5_left') 
left_res = load('~/CS/BAC/CMRNA5/alignTmp_101_CMRNA5_left_finalRes');

% reads
[u_right,n_right]=distReads_CMRNA5(right,right_res);
[u_left,n_left]=distReads_CMRNA5(left,left_res);

figure(1);clf
subplot(2,1,1)
plot(u_right,n_right./sum(n_right),'.','markersize',10)
hold on
plot(u_left,n_left./sum(n_left),'r.','markersize',10)
set(gca,'fontweight','bold','fontsize',15)
title('coverage along the gene')
subplot(2,1,2)
plot(u_right,n_right,'.','markersize',10)
hold on
plot(u_left,n_left,'r.','markersize',10)
set(gca,'fontweight','bold','fontsize',15)
title('coverage along the gene')
print -dpdf ~/CS/BAC/CMRNA5/coverage_along_gene_both_legs


figure(2)
subplot(2,1,1)
plot(u_right,n_right,'.','markersize',10)
set(gca,'ylim',[0 4000])
subplot(2,1,2)
plot(u_left,n_left,'r.','markersize',10)
set(gca,'ylim',[0 4000])


u_all = unique([u_right,u_left]);
n_all = zeros(length(u_all),1);
for i=1:length(n_all)
  n_all(i) = 0;
  a = find(u_right==u_all(i));
  if ~isempty(a)
    n_all(i) = n_all(i)+n_right(a);
  end
  
  a = find(u_left==u_all(i));
  if ~isempty(a)
    n_all(i) = n_all(i)+n_left(a);
  end
end


figure
plot(u_all,n_all,'.')


P_right = round(right_res.POS);
P_left = round(left_res.POS);

sum(right.uniqueReads_length(find(P_right==1248)))
sum(left.uniqueReads_length(find(P_left==1248)))


all = load('~/CS/BAC/CMRNA5/dataCMRNA5') 
all_res = load('~/CS/BAC/CMRNA5/alignTmp_101_CMRNA5_finalRes');


[junk,i1,i2] = intersect(all.uniqueReads,right.uniqueReads,'rows');

all.uniqueReads(i1(1),:)
right.uniqueReads(i2(1),:)

for i=1:length(i1)
  all_res.POS(i1(i))
  right_res.POS(i2(i))
  pause
end


