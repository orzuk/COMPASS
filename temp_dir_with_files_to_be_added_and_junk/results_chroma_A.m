function results_chroma_A(fileName,s1_start,s2_start,n,pos_xl)

load ~/CS/BAC/12Samples/validation/scf_files_drosophila/insilocopcr_sanger_drosophila


basicDirName = '~/CS/BAC/12Samples/validation/scf_files_drosophila/';
ampSequenceName = Sequence_amp_aceto1;
ampHeader = Header_amp_aceto1;

tmp_ind_W_454 = [4452 156437 274161 92777 325754 4421];
tmp_ind_W_Solexa = [65495 274161 4452 58371];


sameNumbers = [4452 4452;274161 274161]; % Solexa,454

list_Solexa_S1_S4 = load('~/CS/BAC/results_for_figure_s1_s4/listOfBAC_Solexa_S1_S4','zz');
list_454_S1_S4 = load('~/CS/BAC/results_for_figure_s1_s4/listOfBAC_454_S1_S4','zz_454');

% use function
results_chroma_in(basicDirName,fileName,ampSequenceName,ampHeader,s1_start,s2_start,n,pos_xl,tmp_ind_W_454,tmp_ind_W_Solexa,list_454_S1_S4,list_Solexa_S1_S4,sameNumbers)

nameForPrint = fileName;
nameForPrint = strrep(nameForPrint,'.scf','')
print('-dpdf',['~/CS/BAC/12Samples/validation/scf_files_drosophila/figs/fig_',nameForPrint])


