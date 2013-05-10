clear

[Header_V2, Sequence_V2] = fastaread('~/CS/BAC/twin/data/V2.fasta');
[Header_V6, Sequence_V6] = fastaread('~/CS/BAC/twin/data/V6.fasta');


l_v2 = zeros(1,length(Sequence_V2));
for i=1:length(Sequence_V2)
  l_v2(i) = length(Sequence_V2{i});
end

l_v6 = zeros(1,length(Sequence_V6));
for i=1:length(Sequence_V6)
  l_v6(i) = length(Sequence_V6{i});
end



v2_tk = find(l_v2>=200);
Sequence_V2_tk = cell(1,length(v2_tk));
for i=1:length(v2_tk)
  Sequence_V2_tk{i} = Sequence_V2{v2_tk(i)}(1:200);
end

v6_tk = find(l_v6>=58);
Sequence_V6_tk = cell(1,length(v6_tk));
for i=1:length(v6_tk)
  Sequence_V6_tk{i} = Sequence_V6{v6_tk(i)}(1:58);
end




% create database of v2 and v6
%1. find all aligned
%2. create 100bp reads
% run original
% run 100bp


% find those with non ACGT

