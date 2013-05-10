function results_chroma_W(fileName,s1_start,s2_start,n,pos_xl)

load ~/CS/BAC/12Samples/validation/scf_files_drosophila/insilocopcr_sanger_drosophila


basicDirName = '~/CS/BAC/12Samples/validation/scf_files_drosophila/';
ampSequenceName = Sequence_amp_wolb2;
ampHeader = Header_amp_wolb2;

tmp_ind_W_454 = [275390 273974 258983 52332 117492 130782 244093 91297 91789 170818];
tmp_ind_W_Solexa = [108447 564629]

sameNumbers = [108447 52332;]; % Solexa,454

list_Solexa_S1_S4 = load('~/CS/BAC/results_for_figure_s1_s4/listOfBAC_Solexa_S1_S4','zz');
list_454_S1_S4 = load('~/CS/BAC/results_for_figure_s1_s4/listOfBAC_454_S1_S4','zz_454');

% use function
results_chroma_in(basicDirName,fileName,ampSequenceName,ampHeader,s1_start,s2_start,n,pos_xl,tmp_ind_W_454,tmp_ind_W_Solexa,list_454_S1_S4,list_Solexa_S1_S4,sameNumbers)

nameForPrint = fileName;
nameForPrint = strrep(nameForPrint,'.scf','')
print('-dpdf',['~/CS/BAC/12Samples/validation/scf_files_drosophila/figs/fig_',nameForPrint])


