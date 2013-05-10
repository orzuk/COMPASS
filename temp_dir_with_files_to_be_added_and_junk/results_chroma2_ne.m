function results_chroma2_ne(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)

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

paperPosition = [0.7    0.01    0.29    0.91];
thresh = 250;
resData = results_chroma2_in(basicDirName,fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,tmp_ind_W_454,tmp_ind_W_Solexa,sameNumbers,paperPosition,thresh);

nameForPrint = fileName;
nameForPrint = strrep(nameForPrint,'.scf','')

%keyboard
print('-dpdf',['~/CS/BAC/12Samples/validation/Human_scf_files/figs/fig2_',nameForPrint])

