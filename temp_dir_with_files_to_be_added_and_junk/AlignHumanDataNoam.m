clear
load ~/CS/mFilesBAC/typicalseqs
auxData = struct;
auxData.brFlag = 0;


if auxData.brFlag
  userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
  disp('run in Br')
else
  userDir = getuserdir;
end

 
auxData.basicSeqKey = [userDir,'/CS/BAC/datNoNonACGT/keyNoNonACGT'];
auxData.basicSeqNameDir = [userDir,'/CS/BAC/datNoNonACGT/packed64/'];

disp('loading basicseqkey');
load(auxData.basicSeqKey,'len_uni');

NUMDBSEQS=410849;
USEREVERSECOMPLEMENT=0;
MAXIDEN=1000;

cppos=zeros(length(typseq),1);
ppos=zeros(MAXIDEN,length(typseq));
numiden=zeros(length(typseq),1);
disp('starting');


part =1:NUMDBSEQS/5:NUMDBSEQS;
part = round(part);
part(end) = NUMDBSEQS+1;

matlabpool open 4

parfor i=1:4
   cppos=zeros(length(typseq),1);
   ppos=zeros(MAXIDEN,length(typseq));
   numiden=zeros(length(typseq),1);
   % each part
   AlignHumanDataPart(i,numiden,cppos,ppos,[part(i),part(i+1)-1],typseq,len_uni,USEREVERSECOMPLEMENT,auxData);
 end
end


for i=1:4
  p{i} = load(['~/CS/BAC/AlignHumanData_output_',num2str(i)]);
end

save ~/CS/BAC/AlignHumanData_compiled p
