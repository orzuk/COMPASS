tic
% clear all
close all

k=2;
[A, C, G, T, ProbA, ProbC, ProbG, ProbT, Comments, PkINdex, Base] = scfread ('5103_L--F.scf');
in_cut=find(smooth(A+T+G+C,100)<10,1,'first');
A=A(1:in_cut)+rand(in_cut,1)*1e-3;
C=C(1:in_cut)+rand(in_cut,1)*1e-3;
G=G(1:in_cut)+rand(in_cut,1)*1e-3;
T=T(1:in_cut)+rand(in_cut,1)*1e-3;
PkINdex=PkINdex(PkINdex<in_cut);
% traceplot(A, C, G, T)
% PkINdex1=repmat(PkINdex+1,1,2*k+1)+repmat([-k:k],length(PkINdex),1);
PkINdex1=cell(length(PkINdex),1);
PkINdex1{1}=[1+k:PkINdex(3)-k];
for i=2:length(PkINdex)-1
    PkINdex1{i}=[PkINdex(i-1)+k:PkINdex(i+1)-k];
end
PkINdex1{end}=[PkINdex(end-1)+k:PkINdex(end)-k];   

% PkINdex1(PkINdex1<1)=1;
% PkINdex1(PkINdex1>length(A))=length(A);
x=zeros(length(PkINdex1),4);
th_par=1;
for i=1:length(PkINdex1)
    
    if i==100
        i
    end
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


[B,index] = sort(x,2,'descend');
th_xmax=10;
th_first2sec=0.05;
for i=1:length(base_opt)
    if x_max(i)>th_xmax
        first2sec = x(i,index(i,2))/x(i,index(i,1));
        if first2sec>th_first2sec
            base_opt(i,index(i,1))=1;
            base_opt(i,index(i,2))=2;
        else
            base_opt(i,index(i,1))=1;
        end
    end
end
    
seq(sum(base_opt,2)==0) = 'N';    
    
    
    
toc