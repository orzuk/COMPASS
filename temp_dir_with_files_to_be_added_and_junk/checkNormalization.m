a = 0.25;
M = [1 1;0 1];
%y = [1-a/3;a/3];
y = [1-2/3*a;2/3*a]
[x]=runOneGroupOf1000ForCompilationFourth(M,y);
b
Mformer = [0.5 1/3;0 1/3];

[x]=runOneGroupOf1000ForCompilationFourth(Mformer,y)
