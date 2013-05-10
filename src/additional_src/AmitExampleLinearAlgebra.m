A=rand(20,4);
A=[A,A(:,1)+2*A(:,2)];
A=A./repmat(sum(A),20,1);%sum of each col is 1
D=rref_mod2(A);
figure; imagesc(D); colorbar;

x=rand(5,1);x=x/sum(x);%the correct vector sum 1
y=A*x;% the measurement vector

x1=lsqr(A,y); % look at the answer you get,
lsqr(A,y)-x % with respect to the correct, note that elements 3,4 are exact but the dependent (col 1,2,5) are wrong

lsqr(A(:,1:4),y)-x(1:4)%if you remove col 5, it doesnt solve the problem
figure; hold on; % why? because the vector y in wrong
plot(y,A(:,1:4)*x(1:4),'.');grid on; hold on;
% once you remove cols the vector y decreas in all elements, to get the
% right y you need to remove the contribution of the 5th col
plot(y-A(:,5)*x(5),A(:,1:4)*x(1:4),'.r');grid on; hold on;
plot(y,y,'k');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load datTmp.mat

A=full(normalizedBac'*normalizedBac);
D=rref_mod2(A);
rank(A)
rank(full(normalizedBac))


x=zeros(500,1);
x(1:200)=rand(200,1);
x=x(randperm(500));
x=x/sum(x);

y=normalizedBac*x;
x1=lsqr(normalizedBac,y,1e-8,1000);

figure; 
plot(x,x1,'.');grid on; hold on;

ind_corsol=abs(x-x1)<1e-6;

%a=diag(D)==1 & sum(abs(D))==1 & sum(abs(D'))==1;%independend independent
%a=diag(D)==1 & [sum(abs(D))==1]' & [sum(abs(D'))==1]';
A(find(abs(D)>10^-5)) = 1;
a = find(sum(A,2)==1 & sum(A,1)'==1 & diag(A)==1);
%a = find(sum(A,2)==1 & sum(A,1)'==1 & diag(A)==1);

plot(x(a),x1(a),'or');grid on; hold on;
plot(x,x,'k')
figure; imagesc(D); colorbar;

N = null(D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

