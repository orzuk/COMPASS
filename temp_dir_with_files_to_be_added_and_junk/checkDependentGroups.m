clear

A = rand(20,6);

A(:,7) = A(:,1)+2*A(:,2);
A(:,8) = A(:,5)+A(:,6);


A = A./repmat(sum(A),size(A,1),1);
A(end+1,:) = 1;
A = A./repmat(sum(A),size(A,1),1);


x = rand(8,1);
x = x./sum(x);
x(end+1) =1;

B=full(A'*A);
D=rref_mod2(B);


