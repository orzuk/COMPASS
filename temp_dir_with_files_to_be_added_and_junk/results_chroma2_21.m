function results_chroma2_21(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)


basicDirName = '~/CS/BAC/12Samples/validation/Human_scf_files/';

tmp_ind_W_454 = [519904 196718 432587];
tmp_ind_W_Solexa = [362429 205914 493003 413657 370504 407929 399892];
 

sameNumbers = [362429 196718;...
               493003 196718;...
               413657  196718;...
               370504   196718;...
               205914  196718;...
               399892  196718;...
              ]; % Solexa,454

% info of those that were found are part of curr_HeaderFull

% use function
%keyboard
paperPosition = [0.7    0.01    0.29    0.91];
thresh = 250;
resData = results_chroma2_in(basicDirName,fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,tmp_ind_W_454,tmp_ind_W_Solexa,sameNumbers,paperPosition,thresh);

nameForPrint = fileName;
nameForPrint = strrep(nameForPrint,'.scf','')

%keyboard
print('-dpdf',['~/CS/BAC/12Samples/validation/Human_scf_files/figs/fig2_',nameForPrint])



