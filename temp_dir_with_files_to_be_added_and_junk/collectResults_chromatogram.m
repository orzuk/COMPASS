clear

% W

clear all;close all
fileName = '5103_W-R.scf';
s1_start = 60; % start of amplified sequence
s2_start = 50; % start of Sanger major allele
n = 600-s2_start;
pos_xl = [100 500];
results_chroma_W(fileName,s1_start,s2_start,n,pos_xl)


clear all;close all
fileName = '5103_W--R.scf';
s1_start = 60; % start of amplified sequence
s2_start = 25; % start of Sanger major allele
n = 720-s2_start;
pos_xl = [100 600];
results_chroma_W(fileName,s1_start,s2_start,n,pos_xl)


clear all;close all
fileName = '5103_W-F.scf';
s1_start = 60; % start of amplified sequence
s2_start = 30; % start of Sanger major allele
n = 620-s2_start;
pos_xl = [100 500];
results_chroma_W(fileName,s1_start,s2_start,n,pos_xl)

clear all;close all
fileName = '5103_W--F.scf';
s1_start = 60; % start of amplified sequence
s2_start = 25; % start of Sanger major allele
n = 500-s2_start;
pos_xl = [100 400];
results_chroma_W(fileName,s1_start,s2_start,n,pos_xl)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% L

clear all;close all
fileName = '5103_L--R.scf';
s1_start = 100; % start of amplified sequence
s2_start = 40; % start of Sanger major allele
n = 680-s2_start;
pos_xl = [100 600];
results_chroma_L(fileName,s1_start,s2_start,n,pos_xl)


clear all;close all
fileName = '5103_L-F.scf';
s1_start = 60; % start of amplified sequence
s2_start = 50; % start of Sanger major allele
n = 200-s2_start;
pos_xl = [100];
results_chroma_L(fileName,s1_start,s2_start,n,pos_xl)

clear all;close all
fileName = '5103_L--F.scf';
s1_start = 60; % start of amplified sequence
s2_start = 40; % start of Sanger major allele
n = 450-s2_start;
pos_xl = [100 400];
results_chroma_L(fileName,s1_start,s2_start,n,pos_xl)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A

clear all;close all
fileName = '5103_A--R.scf';
s1_start = 1; % start of amplified sequence
s2_start = 40; % start of Sanger major allele
n = 550-s2_start;
pos_xl = [100 400];
results_chroma_A(fileName,s1_start,s2_start,n,pos_xl)

clear all;close all
fileName = '5103_A--F.scf';
s1_start = 20; % start of amplified sequence
s2_start = 40; % start of Sanger major allele
n = 500-s2_start;
pos_xl = [100 400];
results_chroma_A(fileName,s1_start,s2_start,n,pos_xl)

%%%%%%%%%%%%%%%%%%%
% Human

% 21
clear all;close all
fileName = '5103_O721-F.scf';
s1_start = 30; % start of amplified sequence
s2_start = 30; % start of Sanger major allele
n = 750-s2_start;
pos_xl = [100 400];
results_chroma_21(fileName,s1_start,s2_start,n,pos_xl)



clear all;close all
fileName = '5103_O721-R.scf';
s1_start = 30; % start of amplified sequence
s2_start = 40; % start of Sanger major allele
n = 800-s2_start;
pos_xl = [100 400];
results_chroma_21(fileName,s1_start,s2_start,n,pos_xl)

clear all;close all
fileName = '5103_O1021-F.scf';
s1_start = 30; % start of amplified sequence
s2_start = 30; % start of Sanger major allele
n = 300-s2_start;
pos_xl = [100 400];
results_chroma_21(fileName,s1_start,s2_start,n,pos_xl)

clear all;close all
fileName = '5103_S721-F.scf';
s1_start = 100; % start of amplified sequence
s2_start = 30; % start of Sanger major allele
n = 800-s2_start;
pos_xl = [100 400];
results_chroma_21(fileName,s1_start,s2_start,n,pos_xl)

clear all;close all
fileName = '5103_S1021-F.scf';
s1_start = 100; % start of amplified sequence
s2_start = 30; % start of Sanger major allele
n = 400-s2_start;
pos_xl = [100 400];
results_chroma_21(fileName,s1_start,s2_start,n,pos_xl)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 51
clear all;close all
fileName = '5103_O751-F.scf';
s1_start = 10; % start of amplified sequence
s2_start = 40; % start of Sanger major allele
n = 400-s2_start;
pos_xl = [100 400];
results_chroma_51(fileName,s1_start,s2_start,n,pos_xl)

clear all;close all
fileName = '5103_O1051-F.scf';
s1_start = 10; % start of amplified sequence
s2_start = 40; % start of Sanger major allele
n = 650-s2_start;
pos_xl = [100 400];
results_chroma_51(fileName,s1_start,s2_start,n,pos_xl)

clear all;close all
fileName = '5103_O1051-R.scf';
s1_start = 50; % start of amplified sequence
s2_start = 50; % start of Sanger major allele
n = 350-s2_start;
pos_xl = [100 400];
results_chroma_51(fileName,s1_start,s2_start,n,pos_xl)

clear all;close all
fileName = '5103_S751-F.scf';
s1_start = 100; % start of amplified sequence
s2_start = 50; % start of Sanger major allele
n = 700-s2_start;
pos_xl = [100 400];
results_chroma_51(fileName,s1_start,s2_start,n,pos_xl)

clear all;close all
fileName = '5103_S751-R.scf';
s1_start = 10; % start of amplified sequence
s2_start = 50; % start of Sanger major allele
n = 800-s2_start;
pos_xl = [100 400];
results_chroma_51(fileName,s1_start,s2_start,n,pos_xl)

clear all;close all
fileName = '5103_S1051-F.scf';
s1_start = 100; % start of amplified sequence
s2_start = 40; % start of Sanger major allele
n = 800-s2_start;
pos_xl = [100 400];
results_chroma_51(fileName,s1_start,s2_start,n,pos_xl)

clear all;close all
fileName = '5103_S1051-R.scf';
s1_start = 100; % start of amplified sequence
s2_start = 50; % start of Sanger major allele
n = 730-s2_start;
pos_xl = [100 400];
results_chroma_51(fileName,s1_start,s2_start,n,pos_xl)

clear all;close all
fileName = '5103_S1051-R.scf';
s1_start = 100; % start of amplified sequence
s2_start = 50; % start of Sanger major allele
n = 730-s2_start;
pos_xl = [100 400];
results_chroma_51(fileName,s1_start,s2_start,n,pos_xl)

% 49
clear all;close all
fileName = '5103_O749-F.scf';
s1_start = 100; % start of amplified sequence
s2_start = 61; % start of Sanger major allele
n = 100-s2_start;
pos_xl = [50];
results_chroma_49(fileName,s1_start,s2_start,n,pos_xl)

clear all;close all
fileName = '5103_S1049-F.scf';
s1_start = 25; % start of amplified sequence
s2_start = 50; % start of Sanger major allele
n = 280-s2_start;
pos_xl = [50];
results_chroma_49(fileName,s1_start,s2_start,n,pos_xl)

% ne
clear all;close all
fileName = '5103_O7NE-F.scf';
s1_start = 10; % start of amplified sequence
s2_start = 50; % start of Sanger major allele
n = 220-s2_start;
pos_xl = [50];
results_chroma_ne(fileName,s1_start,s2_start,n,pos_xl)

clear all;close all
fileName = '5103_O10NE-F.scf';
s1_start = 10; % start of amplified sequence
s2_start = 100; % start of Sanger major allele
n = 400-s2_start;
pos_xl = [50];
results_chroma_ne(fileName,s1_start,s2_start,n,pos_xl)

clear all;close all
fileName = '5103_S7NE-F.scf';
s1_start = 10; % start of amplified sequence
s2_start = 40; % start of Sanger major allele
n = 390-s2_start;
pos_xl = [50];
results_chroma_ne(fileName,s1_start,s2_start,n,pos_xl)

clear all;close all
fileName = '5103_S10NE-F.scf';
s1_start = 10; % start of amplified sequence
s2_start = 40; % start of Sanger major allele
n = 500-s2_start;
pos_xl = [50];
results_chroma_ne(fileName,s1_start,s2_start,n,pos_xl)

454    O7    O10   S7     S10
524892 +na   -      -      -
242939 -na   -      +      +   
365063 +na   -      -      -
217840 +c    -inc   -inc     -inc
569075 +c    -inc   +c      +c

Solexa  O7      O10     S7      S10
569075   +c       +c      -inc   -inc
365063   -na    -na     -na     +na
375840   -c      -inc     +c     -inc
441937   -na    -na     +na     -na
492883   -inc       -inc     -inc     +c

 

O7: 492883 is the same as 217840 over the 450 
O10: 492883 is the same as 217840 over the 450

% does not find the start position. does not find the stop position
