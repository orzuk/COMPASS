function w = l2mu(A,b,numiter)
% w = l2mu(A,b,numiter)
% Given a matrix A and vector b, find w on the simplex which minimizes
% ||Aw-b||_2. Runs multiplicative updates

[m n] = size(A);
eta = 5*sqrt(log(n)/numiter); % Learning rate

w = ones(n,1)/n;
threshold = round(numiter/100); % when to report on progress
fprintf('0   percent done');
for i=1:numiter
    j = randi(m);
    w = w.*(1-eta*(1-(A(j,:)*w-b(j))*A(j,:)*w)); % replace g by 1. Wrong idea!
    w = w/sum(w);
    if (mod(i,threshold)==0) 
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%3d percent done',floor(100*i/numiter))
    end;
end;
fprintf('\n');





