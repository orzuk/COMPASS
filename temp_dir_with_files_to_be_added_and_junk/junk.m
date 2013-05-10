z = redM(:,s{6}); D = rref_mod2(z'*z);figure(1);imagesc(D);figure(2);[x,y] = find(z);imagesc(z(x,:))

zz = full(z(x,:));

k = 8;
rest = 1:size(zz,2);
rest(k) = [];

linsolve(zz(:,rest),zz(:,k))
