load ~/CS/BAC/datTmp

A=full(normalizedBac'*normalizedBac);
D=rref_mod2(A);
rank(A)
rank(full(normalizedBac))


x=zeros(500,1);
x(1:200)=rand(200,1);
x=x(randperm(500));
x=x/sum(x);

y=normalizedBac*x;

[u x_min x_max] = set_solution_bounds(A, x, y, 1)
cvx_begin
 variable s(null_dim)
 minimize(  g*s  );
 subject to
 (-null_basis)*s >= x;
cvx_end
