function results_chroma_L(fileName,s1_start,s2_start,n,pos_xl)

load ~/CS/BAC/12Samples/validation/scf_files_drosophila/insilocopcr_sanger_drosophila


basicDirName = '~/CS/BAC/12Samples/validation/scf_files_drosophila/';
ampSequenceName = Sequence_amp_lacto;
ampHeader = Header_amp_lacto;

tmp_ind_W_454 = [546229 284270 587695 106977 161334 583941 534907 290616 15163 548471 257487 551236];
tmp_ind_W_Solexa = [282170 73879 560153 161334 542366];

sameNumbers = [282170 284270;73879 284270;161334 161334;]; % Solexa,454

list_Solexa_S1_S4 = load('~/CS/BAC/results_for_figure_s1_s4/listOfBAC_Solexa_S1_S4','zz');
list_454_S1_S4 = load('~/CS/BAC/results_for_figure_s1_s4/listOfBAC_454_S1_S4','zz_454');

% use function
results_chroma_in(basicDirName,fileName,ampSequenceName,ampHeader,s1_start,s2_start,n,pos_xl,tmp_ind_W_454,tmp_ind_W_Solexa,list_454_S1_S4,list_Solexa_S1_S4,sameNumbers)

nameForPrint = fileName;
nameForPrint = strrep(nameForPrint,'.scf','')
print('-dpdf',['~/CS/BAC/12Samples/validation/scf_files_drosophila/figs/fig_',nameForPrint])


