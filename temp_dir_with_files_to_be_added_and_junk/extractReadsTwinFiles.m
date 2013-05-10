function [reads,non]=extractReadsTwinFiles(Sequence)
%keyboard
L = zeros(length(Sequence),1);
for i=1:length(Sequence)
  L(i) = length(Sequence{i});
end

I = find(L==101);
reads = Sequence(I);

%keyboard
non = zeros(length(reads),1);
for i=1:length(reads)
  if ~isempty(find(reads{i}~='A' & reads{i}~='C' & reads{i}~='G' & reads{i}~='T'))
    non(i) = 1;
  end
end
non = find(non);
size(non)

