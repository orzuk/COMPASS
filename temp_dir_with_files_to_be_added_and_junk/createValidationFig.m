% W
clear all
primerName = 'wolb2';
load(['~/CS/BAC/12Samples/validation/scf_files_drosophila/resModified_',primerName])

th_sec2first = 0.03;

close all
fileName = '5103_W--F.scf';
AmplifiedSeq_start = 60; % start of amplified sequence
SangerSeq_start = 25; % start of Sanger major allele
lengthValidSanger =480-SangerSeq_start;
pos_xl = [100 400];
%results_chroma4_W(fileName,cData,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,th_sec2first)

results_chroma5_W(fileName,cData,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,th_sec2first)

close all
fileName = '5103_W--R.scf';
AmplifiedSeq_start = 60; % start of amplified sequence
SangerSeq_start = 25; % start of Sanger major allele
lengthValidSanger =680-SangerSeq_start;
pos_xl = [100 600];
results_chroma4_W(fileName,cData,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,th_sec2first)

%%%%%%%%%%%
clear
primerName = 'aceto1';
load(['~/CS/BAC/12Samples/validation/scf_files_drosophila/resModified_',primerName])
%load(['~/CS/BAC/12Samples/validation/scf_files_drosophila/res_',primerName])


th_sec2first = 0.03;  
close all
fileName = '5103_A--R.scf';
AmplifiedSeq_start = 1; % start of amplified sequence
SangerSeq_start = 40; % start of Sanger major allele
lengthValidSanger = 550-SangerSeq_start;
pos_xl = [100 400];
results_chroma4_A(fileName,cData,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,th_sec2first)


close all
filefName = '5103_A--F.scf';
AmplifiedSeq_start = 20; % start of amplified sequence
SangerSeq_start = 40; % start of Sanger major allele
lengthValidSanger = 450-SangerSeq_start;
pos_xl = [100 400];
results_chroma4_A(fileName,cData,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,th_sec2first)







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% W
clear
fileName = '5103_W--F.scf';
nameForPrint = fileName;
nameForPrint = strrep(nameForPrint,'.scf','');
errorF = load(['~/CS/BAC/12Samples/validation/scf_files_drosophila/figs/errorsMinor_',nameForPrint]);
fileName = '5103_W--R.scf';
nameForPrint = fileName;
nameForPrint = strrep(nameForPrint,'.scf','');
errorR = load(['~/CS/BAC/12Samples/validation/scf_files_drosophila/figs/errorsMinor_',nameForPrint]);

[junk,i1,i2] = intersect(errorF.reordered_cData.Numbers,errorR.reordered_cData.Numbers);

E = NaN*ones(length(junk),2);
for i=1:length(junk)
  currSeqF = errorF.errorsMinorOverSequence{i1(i)};
  currSeqR = errorR.errorsMinorOverSequence{i2(i)};
  
  eF = find(currSeqF(1,:)-currSeqF(2,:) & currSeqF(2,:)~='P' &  currSeqF(2,:)=='N');
  eR = find(currSeqR(1,:)-currSeqR(2,:) & currSeqR(2,:)~='P' &  currSeqR(2,:)=='N');

  E(i,:) = [junk(i) length(unique([eF,eR]))];
  length(find(currSeqF(2,:)~='P' |   currSeqR(2,:)~='P'))

end

% length intersection of F and R:

% load the table file and reorder E accordingly

[junk,junk1,raw] = xlsread('~/CS/BAC/12Samples/validation/scf_files_drosophila/figs/WolbachiaTable.xls');

add = zeros(size(E,1),2);
k = 1;
for i=1:size(raw,1)
  if ~isempty(raw{i,1}) && isnumeric(raw{i,1})
    [ind] = find(E(:,1)==raw{i,1});
    add(k,:) = [raw{i,1},E(ind,2)];
    k = k+1;
  end
end


fid = fopen(['~/CS/BAC/12Samples/validation/scf_files_drosophila/figs/W_S2.xls'],'w');
fprintf(fid,'id\t total number of minor errors\n');
for i=1:size(add,1)
  fprintf(fid,'%d\t %d\n',add(i,1),add(i,2));
end
fclose(fid);


%%%%%%%%%%%%%%%%%%5555
% 

%%%%%%%%%%%%%55



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A
clear
fileName = '5103_A--F.scf';
nameForPrint = fileName;
nameForPrint = strrep(nameForPrint,'.scf','');
errorF = load(['~/CS/BAC/12Samples/validation/scf_files_drosophila/figs/errorsMinor_',nameForPrint]);
fileName = '5103_A--R.scf';
nameForPrint = fileName;
nameForPrint = strrep(nameForPrint,'.scf','');
errorR = load(['~/CS/BAC/12Samples/validation/scf_files_drosophila/figs/errorsMinor_',nameForPrint]);

[junk,i1,i2] = intersect(errorF.reordered_cData.Numbers,errorR.reordered_cData.Numbers);

E = NaN*ones(length(junk),2);
for i=1:length(junk)
  currSeqF = errorF.errorsMinorOverSequence{i1(i)};
  currSeqR = errorR.errorsMinorOverSequence{i2(i)};
  
  eF = find(currSeqF(1,:)-currSeqF(2,:) & currSeqF(2,:)~='P' &  currSeqF(2,:)=='N');
  eR = find(currSeqR(1,:)-currSeqR(2,:) & currSeqR(2,:)~='P' &  currSeqR(2,:)=='N');

  E(i,:) = [junk(i) length(unique([eF,eR]))];
  
end

% length intersection of F and R:
length(find(currSeqF(2,:)~='P' |   currSeqR(2,:)~='P'))


% load the table file and reorder E accordingly

[junk,junk1,raw] = xlsread('~/CS/BAC/12Samples/validation/scf_files_drosophila/figs/AcetobacterTable.xls');

add = zeros(size(E,1),2);
k = 1;
for i=1:size(raw,1)
  if ~isempty(raw{i,1}) && isnumeric(raw{i,1})
    [ind] = find(E(:,1)==raw{i,1});
    add(k,:) = [raw{i,1},E(ind,2)];
    k = k+1;
  end
end


fid = fopen(['~/CS/BAC/12Samples/validation/scf_files_drosophila/figs/A_S2.xls'],'w');
fprintf(fid,'id\t total number of minor errors\n');
for i=1:size(add,1)
  fprintf(fid,'%d\t %d\n',add(i,1),add(i,2));
end
fclose(fid);


%%%%%%%%%%%%%%%%%%5555
% 

