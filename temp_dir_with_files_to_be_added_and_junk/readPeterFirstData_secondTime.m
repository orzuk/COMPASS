%%%%%%%%%%%55
% prepare the data
clear
load ~/CS/BAC/full16S/bac16s_full_without_ambiguous

bacseq_len = zeros(1,size(Sequence_uni,2));
for i=1:size(Sequence_uni,2)
  bacseq_len(i) = length(Sequence_uni{i});
end

mx = max(bacseq_len);
Sequence_uni_chararray = char(zeros(mx,size(Sequence_uni,2)));
for i=1:size(Sequence_uni,2)
  Sequence_uni_chararray(1:bacseq_len(i),i) = Sequence_uni{i}';
  Sequence_uni_chararray(bacseq_len(i)+1:end,i) = 'N';
end

save ~/CS/BAC/full16S/bac16s_full_without_ambiguous_charArray Sequence_uni_chararray  bacseq_len 


clear
load ~/CS/BAC/full16S/bac16s_full_without_ambiguous_charArray
x
database1 = Sequence_uni_chararray(:);x
database_size = size(Sequence_uni_chararray);

save ~/CS/BAC/12Samples/Solexa/data/forAlign_full16S/dataForAlign_full16S database1 database_size bacseq_len

%%%%%%%%%%%
%%%%%%%%%%


f = [7000,...
       69592,...
       69611,...
       69620,...
       69656,...
      139471,...
      148441,...
      148442,...
      194145,...
      199987,...
      214360,...
      214361,...
      214362,...
      214557,...
      264886,...
      264890,...
      264893,...
      276733,...
      282902,...
      368938,...
      390968,...
      392371,396770];


seqf = Sequence_uni_chararray(:,f)
Header_f = Header_uni(f);
save ~/CS/BAC/12Samples/validation_Seq seqf Header_f

matrix for all these - how many don't we describe
  add the one which he did 





535   559   559   559   559   557   563   564   531   504   539   539   539   539   560   560   560   558   558   458   561   511   530


% taking only those which are 2*101


clear

[Header, Sequence] = fastaread('~/CS/BAC/Peter/A_Amp_16S.fna');

L = zeros(length(Sequence),1);
for i=1:length(Sequence)
  L(i) = length(Sequence{i});
end

I = find(L==202);
reads = cell(length(I)*2,1);

for i=1:length(I)
  reads{i*2-1,:} = Sequence{I(i)}(1:101);
  reads{i*2,:} = Sequence{I(i)}(102:202);
end

non = zeros(length(reads),1);
for i=1:length(reads)
  if ~isempty(find(reads{i}~='A' & reads{i}~='C' & reads{i}~='G' & reads{i}~='T'))
    non(i) = 1;
  end
end

non = find(non);

reads(non) = [];


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


save ~/CS/BAC/PeterFirstData_secondTime/data/data_PeterFirstData_secondTime uni_reads uniqueReads_length uniqueReads
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% align the reads
clear
%userDir = getuserdir;
userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';

sampleName = 'PeterFirstData_secondTime';

auxData = struct;
auxData.queueName = 'hour';
auxData.readLength = 101;
auxData.searchForwardAndReverseFlag = 0; %genric - pay attanetion if read are from both starnads!!!  
auxData.readsFileName = ['data_PeterFirstData_secondTime'];
auxData.inputDirName = [userDir,'/CS/BAC/PeterFirstData_secondTime/data'];
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
load ~/CS/BAC/Peter/dataPeter
load('~/CS/BAC/Peter/alignTmp_101_Peter_finalRes')

a = find(POS);
length(a)

figure(1)
hist(POS(a),1:1500)
figure(2)

b = unique(POS(a));
num = zeros(length(b),1);
for i=1:length(b)
  c = find(POS(a)==b(i));
  num(i) = sum(uniqueReads_length(a(c)));
end

figure
plot(b,num,'.')

z = round(b);
u = unique(z);
for i=1:length(u)
  n(i) = sum(num(find(z==u(i))));
end

plot(u,n,'.','markersize',10)
set(gca,'fontweight','bold','fontsize',15)
title('coverage along the gene')
print -dpdf ~/CS/BAC/Peter/coverage_along_gene

% plot the distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare the reads - only aligned - run without correction
if 1==2
  uniqueReads = uniqueReads(a,:);
  uniqueReads_length = uniqueReads_length(a);

  save ~/CS/BAC/Peter/uniqueAlignedReadsPeter uniqueReads uniqueReads_length
  clear uniqueReads uniqueReads_length
end

% run unnormalized
runDataCMRNA5_no_eq_no_correction('testCMRNA5_alignedReads_noNor')

% run normalized

% run 90 without normalization

% run 90 with normalization




