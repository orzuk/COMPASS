function results_chroma_49(fileName,s1_start,s2_start,n,pos_xl)

load ~/CS/BAC/12Samples/validation/Human_scf_files/insilocopcr_sanger_human


basicDirName = '~/CS/BAC/12Samples/validation/Human_scf_files/';
ampSequenceName = Sequence_amp_490780;
ampHeader = Header_amp_490780;

tmp_ind_W_454 = [453712 310034 450591 396365 438808 496851 490719 254347 520236 397030 588309 486804 133788 280965 520823 385958];
tmp_ind_W_Solexa = [133788 414748 377679 454080 500036 172279 150299 344193 439046 419776 381508 490780 345714 462062 446695 420637 403439 462166 441866 286967 453712 351798 443886 275755 414931 423251 389536];
 
 



sameNumbers = [254347 490780;...
              ]; % Solexa,454

list_Solexa_Human = load('~/CS/BAC/results_for_figure_o7_s10/listOfBAC_Solexa_Human','zz');
list_454_Human = load('~/CS/BAC/results_for_figure_o7_s10/listOfBAC_454_Human','zz_454');

% use function
results_chroma_in(basicDirName,fileName,ampSequenceName,ampHeader,s1_start,s2_start,n,pos_xl,tmp_ind_W_454,tmp_ind_W_Solexa,list_454_Human,list_Solexa_Human,sameNumbers)

nameForPrint = fileName;
nameForPrint = strrep(nameForPrint,'.scf','')
print('-dpdf',['~/CS/BAC/12Samples/validation/Human_scf_files/figs/fig_',nameForPrint])


