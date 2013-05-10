
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1500
if 1==2
clear all
global Header_750 Sequence_750 header_uni1to350_primers450  seq_uni1to350_primers450 Header_uni Sequence_uni
load ~/CS/BAC/primers750/bac16s_primers750_full_without_ambiguous
load  ~/CS/BAC/primers_startAs750_length350/bac16s_1to350_primers450.mat
load ~/CS/BAC/bacteria_s16_data_uni

createResForSequenceValidation_in('205914')
createResForSequenceValidation_in('490780')
createResForSequenceValidation_in('517821')
createResForSequenceValidation_in('neisseria')
end
%%%%%%%%%%%%
% 21

if 1==2
clear
primerName = '205914';
load(['~/CS/BAC/12Samples/validation/Human_scf_files/res_',primerName])

% do 21
fileName = '5103_O721-F.scf';
SangerSeq_start = 30; % start of Sanger major allele
AmplifiedSeq_start = 1;
lengthValidSanger = 750-SangerSeq_start; % 750 is the length of the Sanger
pos_xl = [100 700]
results_chroma2_21(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)

fileName = '5103_O721-R.scf';
SangerSeq_start = 40; % start of Sanger major allele
AmplifiedSeq_start = 1;
lengthValidSanger = 800-SangerSeq_start; % 750 is the length of the Sanger
pos_xl = [100 700]
results_chroma2_21(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)

fileName = '5103_O1021-F.scf';
SangerSeq_start = 30; % start of Sanger major allele
AmplifiedSeq_start = 1;
lengthValidSanger = 300-SangerSeq_start; % 750 is the length of the Sanger
pos_xl = [100 200]
results_chroma2_21(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)

fileName = '5103_S721-F.scf';
AmplifiedSeq_start = 100;
SangerSeq_start = 30; % start of Sanger major allele
lengthValidSanger = 800-SangerSeq_start; % 750 is the length of the Sanger
pos_xl = [100 400];
results_chroma2_21(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)
% did not work

fileName = '5103_S1021-F.scf';
AmplifiedSeq_start = 100;
SangerSeq_start = 30; % start of Sanger major allele
lengthValidSanger = 400-SangerSeq_start; % 750 is the length of the Sanger
pos_xl = [100 400];
results_chroma2_21(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)
% did not work


end

%%%%%%%%%%%
% 51
if 1==2

  clear
primerName = '517821';
load(['~/CS/BAC/12Samples/validation/Human_scf_files/res_',primerName])

close all
fileName = '5103_O751-F.scf';
AmplifiedSeq_start = 10; % start of amplified sequence
SangerSeq_start = 40; % start of Sanger major allele
lengthValidSanger = 400-SangerSeq_start;
pos_xl = [100 400];
results_chroma2_51(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)

close all
fileName = '5103_O1051-F.scf';
AmplifiedSeq_start = 10; % start of amplified sequence
SangerSeq_start = 40; % start of Sanger major allele
lengthValidSanger = 650-SangerSeq_start;
pos_xl = [100 400];
results_chroma2_51(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)


close all
fileName = '5103_O1051-R.scf';
AmplifiedSeq_start = 50; % start of amplified sequence
SangerSeq_start = 50; % start of Sanger major allele
lengthValidSanger = 350-SangerSeq_start;
pos_xl = [100 400];
results_chroma2_51(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)


close all
fileName = '5103_S751-F.scf';
AmplifiedSeq_start = 100; % start of amplified sequence
SangerSeq_start = 50; % start of Sanger major allele
lengthValidSanger = 700-SangerSeq_start;
pos_xl = [100 400];
results_chroma2_51(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)


close all
fileName = '5103_S751-R.scf';
AmplifiedSeq_start = 10; % start of amplified sequence
SangerSeq_start = 50; % start of Sanger major allele
lengthValidSanger = 800-SangerSeq_start;
pos_xl = [100 400];
results_chroma2_51(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)


close all
fileName = '5103_S1051-F.scf';
AmplifiedSeq_start = 100; % start of amplified sequence
SangerSeq_start = 40; % start of Sanger major allele
lengthValidSanger = 800-SangerSeq_start;
pos_xl = [100 400];
results_chroma2_51(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)


close all
fileName = '5103_S1051-R.scf';
AmplifiedSeq_start = 100; % start of amplified sequence
SangerSeq_start = 50; % start of Sanger major allele
lengthValidSanger = 730-SangerSeq_start;
pos_xl = [100 400];
results_chroma2_51(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)


close all
fileName = '5103_S1051-R.scf';
AmplifiedSeq_start = 100; % start of amplified sequence
SangerSeq_start = 50; % start of Sanger major allele
lengthValidSanger = 730-SangerSeq_start;
pos_xl = [100 400];
results_chroma2_51(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)


  
end
%%%%%%%%5 % end 51


% NE
if 1==1

clear
primerName = 'neisseria';
load(['~/CS/BAC/12Samples/validation/Human_scf_files/res_',primerName])

  close all
fileName = '5103_O7NE-F.scf';
AmplifiedSeq_start = 10; % start of amplified sequence
SangerSeq_start = 50; % start of Sanger major allele
lengthValidSanger = 220-SangerSeq_start;
pos_xl = [50];
results_chroma2_ne(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)


close all
fileName = '5103_O10NE-F.scf';
AmplifiedSeq_start = 10; % start of amplified sequence
SangerSeq_start = 100; % start of Sanger major allele
lengthValidSanger = 400-SangerSeq_start;
pos_xl = [50];
results_chroma2_ne(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)

close all
fileName = '5103_S7NE-F.scf';
AmplifiedSeq_start = 10; % start of amplified sequence
SangerSeq_start = 40; % start of Sanger major allele
lengthValidSanger = 390-SangerSeq_start;
pos_xl = [50];
results_chroma2_ne(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)

close all
fileName = '5103_S10NE-F.scf';
AmplifiedSeq_start = 10; % start of amplified sequence
SangerSeq_start = 40; % start of Sanger major allele
lengthValidSanger = 500-SangerSeq_start;
pos_xl = [50];
results_chroma2_ne(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)

end


