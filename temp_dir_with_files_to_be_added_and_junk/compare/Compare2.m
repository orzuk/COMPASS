% compare original solution and reconstruction
% by joining sequences into groups so that min(k-mer distance(seq from orig
% group, seq from reconstructed group))<dist
%
% input:
% mat - a matrix containing >0 if k-mer i is in sequence j
% set1 - the columns of mat which contain the sequences in the original
% mixture (length n1)
% freq1 - (length n1) the frequency of each sequence present in the mixture
% set2 - the columns of mat which contain the sequences in the
% reconstructed mixture (length n2)
% freq2 - (length n2) the frequency of each sequence present in the
% reconstruction
% dist - the threshold distance (in #of k-mers) under which 2 sequences are
% defined similar (51 for 1 mismatch in 50-mers)
%
% output:
% freqseq: a 2*N matrix of group frequencies
% row1 - the frequency of the group in the original mixture
% row2 - the frequency of the group in the reconstructed mixture
function [freqset]=Compare2(mat,set1,freq1,set2,freq2,dist)

similar=zeros(length(set2),1);
used=zeros(length(set2),1);
freqset=[];
for a=1:length(set1)
    totsimilar=0;
    for b=1:length(set2)
        numsim=length(find(mat(:,set1(a)).*mat(:,set2(b))>0));
        if (min(nnz(mat(:,set1(a))),nnz(mat(:,set2(b))))-numsim<dist)
            similar(totsimilar+1)=b;
            totsimilar=totsimilar+1;
        end
    end
    if (totsimilar>0)
        usedorig=[];
        addfreq=0;
        for c=1:totsimilar
            if (used(similar(c))>0)
                usedorig=[usedorig used(similar(c))];
            else
                addfreq=addfreq+freq2(similar(c));
                used(similar(c))=a;
            end
        end
        if (isempty(usedorig))
            freqset=[freqset;freq1(a),addfreq];
        else
%            disp(['not empty ' num2str(a)]);
            freqset=[freqset;0,0];
            usedorig=unique(usedorig);
            minjoin=min(usedorig);
            usedorig(find(usedorig==minjoin))=[];
            used(find(used==a))=minjoin;
            freqset(minjoin,1)=freqset(minjoin,1)+freq1(a);
            freqset(minjoin,2)=freqset(minjoin,2)+addfreq;
            for c=1:length(usedorig)
                disp(['Joining orig sets ' num2str(minjoin) ' and ' num2str(usedorig(c))]);
                freqset(minjoin,:)=freqset(minjoin,:)+freqset(usedorig(c),:);
                freqset(usedorig(c),1)=0;
                freqset(usedorig(c),2)=0;
                used(find(used==usedorig(c)))=minjoin;
            end
        end
    else
        freqset=[freqset;freq1(a),0];
        disp(['No similar to ' num2str(a)]);
    end
end
notused=find(used==0);
for a=1:length(notused)
    freqset=[freqset;0,freq2(notused(a))];
end
