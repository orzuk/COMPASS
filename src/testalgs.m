% Script to test optimization algorithms
% Creates a "random" matrix A, where each column is a sparse distribution
% vector over possible readings. Optimization is then ran with respect to
% the "true" frequencies vector wtrue, which is set to be wtrue(2:5)=0.25, 
% and 0 otherwise.

%%%% PARAMETERS %%%%%%%
n = 3000; % Number of species and possible readings
numseq = 20; % Number of sequences per species
numiter = 30000; % Number of iterations

%%%%%% Creating Matrix A %%%%%%
i = zeros(numseq*n,1);
j = zeros(numseq*n,1);
s = zeros(numseq*n,1);
counter = 0;
fprintf('creating matrix A...\n');
fprintf('0   percent done');
threshold = round(n/100); % when to report on progress
for c=1:n
    j(counter+1:counter+numseq) = c;
    a1 = randperm(n);
    i(counter+1:counter+numseq) = a1(1:numseq); % Choose non-zero entries
    a2 = randn(numseq,1); a2=abs(a2); a2=a2/sum(a2); % Choose random distribution
    s(counter+1:counter+numseq) = a2;
    counter = counter + numseq;
    if (mod(c,threshold)==0) 
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%3d percent done',floor(100*c/n))
    end;    
end;
clear a1 a2;
A = sparse(i,j,s,n,n);
clear i j s
fprintf('\nfinished creating A.\n');

%%%%%% Creating Coordinate Descent %%%%%%
fprintf('Running coordinate descent algorithm...\n');
wtrue = zeros(n,1); wtrue(2:5) = 0.25;
clear dist;
w = l2cd(A,A*wtrue,numiter);
fprintf('values of w on 4 real species: %f %f %f %f\n',w(2),w(3),w(4),w(5));
fprintf('maximal value of w on other species: %f\n',max(w(1),max(w(6:end))));
%fprintf('Distance from truth in L1 norm: %f\n',sum(abs(w-wtrue)));
%fprintf('Distance from truth in L2 norm: %f\n',norm(w-wtrue));

%%%%%% Creating Multiplicative Updates %%%%%%
fprintf('Running multiplicative updates algorithm...\n');
wtrue = zeros(n,1); wtrue(2:5) = 0.25;
clear dist;
w = l2mu(A,A*wtrue,numiter); % This should change to multiplicative updates !!! 
fprintf('values of w on 4 real species: %f %f %f %f\n',w(2),w(3),w(4),w(5));
fprintf('maximal value of w on other species: %f\n',max(w(1),max(w(6:end))));
%fprintf('Distance from truth in L1 norm: %f\n',sum(abs(w-wtrue)));
%fprintf('Distance from truth in L2 norm: %f\n',norm(w-wtrue));