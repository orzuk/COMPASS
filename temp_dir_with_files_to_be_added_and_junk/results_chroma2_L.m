function results_chroma2_L(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)

load ~/CS/BAC/12Samples/validation/scf_files_drosophila/insilocopcr_sanger_drosophila


basicDirName = '~/CS/BAC/12Samples/validation/scf_files_drosophila/';
ampSequenceName = Sequence_amp_lacto;
ampHeader = Header_amp_lacto;

tmp_ind_W_454 = [546229 284270 587695 106977 161334 583941 534907 290616 15163 548471 257487 551236];
tmp_ind_W_Solexa = [282170 73879 560153 161334 542366];

sameNumbers = [282170 284270;73879 284270;161334 161334;]; % Solexa,454

paperPosition = [0.8    0.001    0.19    0.96];
fs = 6;
thresh = 380;
disp('v3');pause
resData = results_chroma3_in(basicDirName,fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,tmp_ind_W_454,tmp_ind_W_Solexa,sameNumbers,paperPosition,thresh,fs);

nameForPrint = fileName;
nameForPrint = strrep(nameForPrint,'.scf','')

%keyboard
print('-dpdf',['~/CS/BAC/12Samples/validation/scf_files_drosophila/figs/fig2_',nameForPrint])
