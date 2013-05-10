function results_chroma2_W(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)

load ~/CS/BAC/12Samples/validation/scf_files_drosophila/insilocopcr_sanger_drosophila


basicDirName = '~/CS/BAC/12Samples/validation/scf_files_drosophila/';
ampSequenceName = Sequence_amp_wolb2;
ampHeader = Header_amp_wolb2;

tmp_ind_W_454 = [275390 273974 258983 52332 117492 130782 244093 91297 91789 170818];
tmp_ind_W_Solexa = [108447 564629]

sameNumbers = [108447 52332;]; % Solexa,454

paperPosition = [0.75    0    0.24    0.97];
fs = 6;
thresh = 400;
resData = results_chroma2_in(basicDirName,fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,tmp_ind_W_454,tmp_ind_W_Solexa,sameNumbers,paperPosition,thresh,fs);

nameForPrint = fileName;
nameForPrint = strrep(nameForPrint,'.scf','')

%keyboard
print('-dpdf',['~/CS/BAC/12Samples/validation/scf_files_drosophila/figs/fig2_',nameForPrint])

