b = ind_bac_in_mix;


basicSeqNameDirOld = ['~/CS/BAC/datNoNonACGT/unpacked/'];
tic
[normalizedBac1 values1] = prepareGroupOf1000DistributedSequenceFiles2(readLength,b,basicSeqNameDirOld,basicSeqKey); 
toc

tic
[normalizedBac values] = prepareGroupOf1000DistributedSequenceFilesOrFourth(readLength,b,basicSeqNameDir,basicSeqKey); 
toc

unp_values1 = nt2int(char(values1));
unp_values = unpack_seqs(values,readLength,64);

[junk,z1,z2] = intersect(unp_values,unp_values1,'rows');


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

