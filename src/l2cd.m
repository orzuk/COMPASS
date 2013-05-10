function w = l2cd(A,b,numiter)
% w = l2cd(A,b,numiter)
% Given a matrix A and vector b, find w on the simplex which minimizes
% ||Aw-b||_2. Runs stochastic coordinate descent for numiter iterations

[m n] = size(A);
w = zeros(n,1); w(randi(n))=1;
threshold = round(numiter/100); % when to report on progress
fprintf('0   percent done');
for i=1:numiter
    j = randi(n);
    p = A*w;
    r1 = p-b;
    r2 = p-A(:,j);
    s = (r1'*r2)/(r2'*r2);
    s = max(0,min(1,s));
    w = (1-s)*w; w(j)=w(j)+s;
    if (mod(i,threshold)==0) 
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%3d percent done',floor(100*i/numiter))
    end;
end;
fprintf('\n');





