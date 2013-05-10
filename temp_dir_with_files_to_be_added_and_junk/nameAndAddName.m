function [w,addName]=nameAndAddName(dirName,nameToTest)

a = findstr(nameToTest,'bac_dist_flag');
tmpName = nameToTest(a:end);
a = find(tmpName=='_');
bdf = tmpName(a(3)+1:a(4)-1);


a = findstr(nameToTest,'Nbac_in_mixture');
tmpName = nameToTest(a:end);
a = find(tmpName=='_');
nb = tmpName(a(3)+1:a(4)-1);

a = findstr(nameToTest,'npower');
tmpName = nameToTest(a:end);
a = find(tmpName=='_');
np = tmpName(a(1)+1:a(2)-1);
if np=='5'
  np = '05';
end

a = findstr(nameToTest,'readlen');
tmpName = nameToTest(a:end);
a = find(tmpName=='_');
readlen = tmpName(a(1)+1:a(2)-1);

a = findstr(nameToTest,'Nreads');
tmpName = nameToTest(a:end);
a = find(tmpName=='_');
numReads = tmpName(a(1)+1:a(2)-1);

addName = ['_Nbacmix_',nb,'_Nread_',numReads,'_Readlen_',readlen,'_Npower_',np,'_bacdistflag_',bdf,'_NoReads'];

w = ['scp ',...
     ' /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/',dirName,'/structFor_',nameToTest,'.mat',...
     ' /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/',dirName,'/',nameToTest,addName,'.mat'...
     ' shental@sol.cslab.openu.ac.il:~/tmp/'];
