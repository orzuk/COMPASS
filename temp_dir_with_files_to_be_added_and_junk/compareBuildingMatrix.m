
clear
%userDir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen';
userDir = getuserdir;

addpath([userDir,'/CS/mFilesBAC'])


basicSeqKey = [userDir,'/CS/BAC/dat450000/key450000'];
basicSeqNameDir = [userDir,'/CS/BAC/dat450000/'];
readLength = 50;

%load ~/CS/BAC/listOfDoubles
%b = 1:800;
%b = setdiff(b,d);

c = load('~/b') 
b = c.b(49);


b=1:2000;
tic
[normalizedBac values1] = prepareGroupOf1000DistributedSequenceFiles2(readLength,b,basicSeqNameDir,basicSeqKey); 
toc

%load([userDir,'/CS/BAC/s16_data_uni_packed64'])
%load('~/CS/BAC/s16_data_uni_packed')
Sequence_packed1 = Sequence_packed64(b);

tic
[ReadBySpeciesMat values] = BuildMixingMatrixFromSequences(readLength, Sequence_packed1, len_uni(b)); % find for last 1000
toc
size(normalizedBac)
size(ReadBySpeciesMat)


figure(1)
imagesc(values1)

figure(2)
m = abs(int2nt(unpack_seqs(values,50,64)));
m = uint8(m);
imagesc(m)

find(values1-m(:,[32:-1:1,33:end]))

c = sortrows(values1)-sortrows(m);


int2nt(unpack_seqs(Sequence_packed1{1},len_uni(b),64))

for i=b
  s{i} = unpack_seqs(Sequence_packed1{i},len_uni(i),64);
end

a = zeros(size(normalizedBac));
a(find(normalizedBac)) =1;
z = sum(a,2);

a = zeros(size(ReadBySpeciesMat));
a(find(ReadBySpeciesMat)) =1;
w = sum(a,2);

q = unpack_seqs(values(50,:),readLength);
q = int2nt(q);
[junk,i1,i2]=intersect(char(values1),q,'rows');

