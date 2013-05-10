function results_chroma2_51(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)

load ~/CS/BAC/12Samples/validation/Human_scf_files/insilocopcr_sanger_human

basicDirName = '~/CS/BAC/12Samples/validation/Human_scf_files/';
ampSequenceName = Sequence_amp_517821;
ampHeader = Header_amp_517821;

tmp_ind_W_454 = [524641];
tmp_ind_W_Solexa = [517821];
 

sameNumbers = [517821 524641;...
              ]; % Solexa,454

paperPosition = [0.5    0.01    0.49    0.91];
thresh = 150;
resData = results_chroma2_in(basicDirName,fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,tmp_ind_W_454,tmp_ind_W_Solexa,sameNumbers,paperPosition,thresh);

nameForPrint = fileName;
nameForPrint = strrep(nameForPrint,'.scf','')

%keyboard
print('-dpdf',['~/CS/BAC/12Samples/validation/Human_scf_files/figs/fig2_',nameForPrint])


