function results_chroma_21(fileName,s1_start,s2_start,n,pos_xl)

load ~/CS/BAC/12Samples/validation/Human_scf_files/insilocopcr_sanger_human


basicDirName = '~/CS/BAC/12Samples/validation/Human_scf_files/';
ampSequenceName = Sequence_amp_205914;
ampHeader = Header_amp_205914;

tmp_ind_W_454 = [519904 196718 432587];
tmp_ind_W_Solexa = [362429 205914 493003 413657 370504 407929 399892];
 

sameNumbers = [362429 196718;...
               493003 196718;...
               413657  196718;...
               370504   196718;...
               205914  196718;...
               399892  196718;...
              ]; % Solexa,454

list_Solexa_Human = load('~/CS/BAC/results_for_figure_o7_s10/listOfBAC_Solexa_Human','zz');
list_454_Human = load('~/CS/BAC/results_for_figure_o7_s10/listOfBAC_454_Human','zz_454');

% use function
results_chroma_in(basicDirName,fileName,ampSequenceName,ampHeader,s1_start,s2_start,n,pos_xl,tmp_ind_W_454,tmp_ind_W_Solexa,list_454_Human,list_Solexa_Human,sameNumbers)

nameForPrint = fileName;
nameForPrint = strrep(nameForPrint,'.scf','')
print('-dpdf',['~/CS/BAC/12Samples/validation/Human_scf_files/figs/fig_',nameForPrint])


