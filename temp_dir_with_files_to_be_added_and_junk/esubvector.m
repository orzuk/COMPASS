function [aa,maxv]=esubvector(a,b,startpos,endpos)
% find the EXACT occurance of a in b (and the match length)
ll=length(a);
score=zeros(ll,1);
for aa=startpos: endpos
    score(aa)=sum(b(aa:aa+ll-1)==a);
end
[maxv aa]=max(score);
