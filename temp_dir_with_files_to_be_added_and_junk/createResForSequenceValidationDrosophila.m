if 1==2
clear all
global Header_750 Sequence_750 header_uni1to350_primers450  seq_uni1to350_primers450 Header_uni Sequence_uni
load ~/CS/BAC/primers750/bac16s_primers750_full_without_ambiguous
load  ~/CS/BAC/primers_startAs750_length350/bac16s_1to350_primers450.mat
load ~/CS/BAC/bacteria_s16_data_uni

createResForSequenceValidationDrosophila_in('wolb2');
createResForSequenceValidationDrosophila_in('aceto1');
createResForSequenceValidationDrosophila_in('lacto');
 
end

if 1==1
% W

clear
primerName = 'wolb2';
load(['~/CS/BAC/12Samples/validation/scf_files_drosophila/res_',primerName])

th_sec2first = 0.02;

close all
fileName = '5103_W-R.scf';
AmplifiedSeq_start = 60; % start of amplified sequence
SangerSeq_start = 50; % start of Sanger major allele
lengthValidSanger =600-SangerSeq_start;
pos_xl = [100 500];
results_chroma3_W(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,th_sec2first)


close all
fileName = '5103_W--R.scf';
AmplifiedSeq_start = 60; % start of amplified sequence
SangerSeq_start = 25; % start of Sanger major allele
lengthValidSanger =680-SangerSeq_start;
pos_xl = [100 600];
results_chroma3_W(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,th_sec2first)



close all
fileName = '5103_W-F.scf';
AmplifiedSeq_start = 60; % start of amplified sequence
SangerSeq_start = 30; % start of Sanger major allele
lengthValidSanger =620-SangerSeq_start;
pos_xl = [100 500];
results_chroma3_W(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,th_sec2first)

close all
fileName = '5103_W--F.scf';
AmplifiedSeq_start = 60; % start of amplified sequence
SangerSeq_start = 25; % start of Sanger major allele
lengthValidSanger =480-SangerSeq_start;
pos_xl = [100 400];
results_chroma3_W(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,th_sec2first)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end

% lacto
if 1==1
clear
primerName = 'lacto';
load(['~/CS/BAC/12Samples/validation/scf_files_drosophila/res_',primerName])

th_sec2first = 0.02;
close all
fileName = '5103_L--R.scf';
AmplifiedSeq_start = 100; % start of amplified sequence
SangerSeq_start = 40; % start of Sanger major allele
lengthValidSanger = 680-SangerSeq_start;
pos_xl = [100 600];
results_chroma3_L(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,th_sec2first)



close all
fileName = '5103_L-F.scf';
AmplifiedSeq_start = 60; % start of amplified sequence
SangerSeq_start = 50; % start of Sanger major allele
lengthValidSanger = 200-SangerSeq_start;
pos_xl = [100];
results_chroma3_L(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,th_sec2first)

close all
fileName = '5103_L--F.scf';
AmplifiedSeq_start = 60; % start of amplified sequence
SangerSeq_start = 40; % start of Sanger major allele
lengthValidSanger = 450-SangerSeq_start;
pos_xl = [100 400];
results_chroma3_L(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,th_sec2first)


end

% aceto
if 1==1
clear
primerName = 'aceto1';
load(['~/CS/BAC/12Samples/validation/scf_files_drosophila/res_',primerName])

th_sec2first = 0.02;  
close all
fileName = '5103_A--R.scf';
AmplifiedSeq_start = 1; % start of amplified sequence
SangerSeq_start = 40; % start of Sanger major allele
lengthValidSanger = 550-SangerSeq_start;
pos_xl = [100 400];
results_chroma3_A(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,th_sec2first)


close all
fileName = '5103_A--F.scf';
AmplifiedSeq_start = 20; % start of amplified sequence
SangerSeq_start = 40; % start of Sanger major allele
lengthValidSanger = 450-SangerSeq_start;
pos_xl = [100 400];
results_chroma3_A(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl,th_sec2first)



  
end