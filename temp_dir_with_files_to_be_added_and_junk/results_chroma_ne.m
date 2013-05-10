function results_chroma_ne(fileName,s1_start,s2_start,n,pos_xl)

load ~/CS/BAC/12Samples/validation/Human_scf_files/insilocopcr_sanger_human


basicDirName = '~/CS/BAC/12Samples/validation/Human_scf_files/';
ampSequenceName = Sequence_amp_neisseria;
ampHeader = Header_amp_neisseria;

tmp_ind_W_454 = [524892 242939 365063 217840 569075 ];
tmp_ind_W_Solexa = [492883 441937 375840 365063 569075 ];


sameNumbers = [569075 569075;...
              365063 365063;...
              441937 524892;...
              492883 217840]; % Solexa,454

list_Solexa_Human = load('~/CS/BAC/results_for_figure_o7_s10/listOfBAC_Solexa_Human','zz');
list_454_Human = load('~/CS/BAC/results_for_figure_o7_s10/listOfBAC_454_Human','zz_454');

% use function
results_chroma_in(basicDirName,fileName,ampSequenceName,ampHeader,s1_start,s2_start,n,pos_xl,tmp_ind_W_454,tmp_ind_W_Solexa,list_454_Human,list_Solexa_Human,sameNumbers)

nameForPrint = fileName;
nameForPrint = strrep(nameForPrint,'.scf','')
print('-dpdf',['~/CS/BAC/12Samples/validation/Human_scf_files/figs/fig_',nameForPrint])


