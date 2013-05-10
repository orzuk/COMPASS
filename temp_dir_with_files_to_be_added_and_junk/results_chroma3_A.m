function results_chroma3_A(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,th_sec2first)

load ~/CS/BAC/12Samples/validation/scf_files_drosophila/insilocopcr_sanger_drosophila


basicDirName = '~/CS/BAC/12Samples/validation/scf_files_drosophila/';
ampSequenceName = Sequence_amp_aceto1;
ampHeader = Header_amp_aceto1;

tmp_ind_W_454 = [4452 156437 274161 92777 325754 4421];
tmp_ind_W_Solexa = [65495 274161 4452 58371];


sameNumbers = [4452 4452;274161 274161]; % Solexa,454
paperPosition = [0.9    0.001    0.09    0.93];
fs = 10;
thresh = 600;
disp('v3');
results_chroma3_in(basicDirName,fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,tmp_ind_W_454,tmp_ind_W_Solexa,sameNumbers,paperPosition,thresh,fs,th_sec2first);

nameForPrint = fileName;
nameForPrint = strrep(nameForPrint,'.scf','')

%keyboard
print('-dpdf',['~/CS/BAC/12Samples/validation/scf_files_drosophila/figs/fig2_',nameForPrint])


