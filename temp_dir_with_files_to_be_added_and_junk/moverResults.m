function moverResults(outFileDirName)

unix(['cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/;tar cf ',outFileDirName,'_res.tar ',outFileDirName,'/']);
unix(['cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/;tar cf ',outFileDirName,'_d.tar ',outFileDirName,'/']);
unix(['cd /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/;tar cf ',outFileDirName,'_t.tar ',outFileDirName,'/']);

['scp /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/',outFileDirName,'_res.tar /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/dataForSim/',outFileDirName,'_d.tar /seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/tmp/',outFileDirName,'_t.tar shental@sol.cslab.openu.ac.il:~/CS/BAC ']