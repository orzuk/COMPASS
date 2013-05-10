function [seq,base_opt]=find_possible_seq_from_chromatogram_noam3(fileName,th_sec2first)

%keyboard

[A, C, G, T, ProbA, ProbC, ProbG, ProbT, Comments, PkINdex, Base] = scfread(fileName);
%keyboard
minVal = max(min(smooth(A+T+G+C,100))+1,10);
in_cut=find(smooth(A+T+G+C,100)<minVal,1,'first');
q = rand(in_cut,1);
A=A(1:in_cut)+q*1e-3;
C=C(1:in_cut)+q*1e-3;
G=G(1:in_cut)+q*1e-3;
T=T(1:in_cut)+q*1e-3;

PkINdex=PkINdex(PkINdex<in_cut);

k=3; %k=2
PkINdex1=cell(length(PkINdex),1);
PkINdex1{1}=[1+k:PkINdex(3)-k];
for i=2:length(PkINdex)-1
    PkINdex1{i}=[PkINdex(i-1)+k:PkINdex(i+1)-k];
end
PkINdex1{end}=[PkINdex(end-1)+k:PkINdex(end)-k];   

%keyboard
x=zeros(length(PkINdex1),4);
th_par=1;
for i=1:length(PkINdex1)
  
    inds=PkINdex1{i};
    if length(inds)>2
        
        t1=A(inds);
        [pksA] = max(findpeaks(t1,'threshold',th_par));
        [pksA0] = max(findpeaks(t1,'threshold',0));
        
        t1=C(inds);
        [pksC] = max(findpeaks(t1,'threshold',th_par));
        [pksC0] = max(findpeaks(t1,'threshold',0));
        
        t1=G(inds);
        [pksG] = max(findpeaks(t1,'threshold',th_par));
        [pksG0] = max(findpeaks(t1,'threshold',0));
        
        t1=T(inds);
        [pksT] = max(findpeaks(t1,'threshold',th_par));
        [pksT0] = max(findpeaks(t1,'threshold',0));
        
        if isempty(pksA0); pksA0=0; end
        if isempty(pksA); pksA=0; end
        if isempty(pksC0); pksC0=0; end
        if isempty(pksC); pksC=0; end
        if isempty(pksG0); pksG0=0; end
        if isempty(pksG); pksG=0; end
        if isempty(pksT0); pksT0=0; end
        if isempty(pksT); pksT=0; end
        
        % change 9.9.12
        pksA0 = max([pksA0,pksA]);pksA = pksA0;
        pksC0 = max([pksC0,pksC]);pksC = pksC0;
        pksG0 = max([pksG0,pksG]);pksG = pksG0; 
        pksT0 = max([pksT0,pksT]);pksT = pksT0;
        % end change
        
        
        pkmax=max([pksA0,pksC0,pksT0,pksG0]);
        if ~isempty(pksA0)
            if pksA0==pkmax
                x(i,1)=max(pksA0);
            else
                x(i,1)=max(pksA);
            end
        end
        if ~isempty(pksC0)
            if pksC0==pkmax
                x(i,2)=max(pksC0);
            else
                x(i,2)=max(pksC);
            end
        end
        if ~isempty(pksG0)
            if pksG0==pkmax
                x(i,3)=max(pksG0);
            else
                x(i,3)=max(pksG);
            end
        end
        if ~isempty(pksT0)
            if pksT0==pkmax
                x(i,4)=max(pksT0);
            else
                x(i,4)=max(pksT);
            end
        end
    end

end

         
% x = [max(A(PkINdex1),[],2), max(C(PkINdex1),[],2), max(G(PkINdex1),[],2), max(T(PkINdex1),[],2)];
[x_max,im]=max(x,[],2);
xmaxratio=x./ repmat(max(x,[],2),1,4);
base_opt = zeros(size(x));
% base_opt((im -1)*length(PkINdex) + [1:length(PkINdex)]') = true;
tmp = sort(xmaxratio');
tmp = tmp(3,:);
% im_second = 
bases='ACGT';
seq=bases(im);
seq(sum(x,2)==0) = 'N';

%keyboard
[B,index] = sort(x,2,'descend');
th_xmax=10;
%keyboard
for i=1:length(base_opt)
    if x_max(i)>th_xmax
        sec2first = x(i,index(i,2))/x(i,index(i,1));
        if sec2first>th_sec2first
            base_opt(i,index(i,1))=1;
            base_opt(i,index(i,2))=2;
        else
            base_opt(i,index(i,1))=1;
        end
        
        third2first = x(i,index(i,3))/x(i,index(i,1));
        if third2first>th_sec2first
            base_opt(i,index(i,3))=3;
        end
        
        fourth2first = x(i,index(i,4))/x(i,index(i,1));
        if fourth2first>th_sec2first
            base_opt(i,index(i,4))=4;
        end
        
        
        
        
    end
end
    
seq(sum(base_opt,2)==0) = 'N';    

