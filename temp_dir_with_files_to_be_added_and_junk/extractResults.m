function extractResults(outFileDirName)

unix(['cd ~/CS/BAC/;tar xf ',outFileDirName,'_res.tar ']);
unix(['cd ~/CS/BAC/dataForSim/;mv ~/CS/BAC/',outFileDirName,'_d.tar .;',' tar xf ',outFileDirName,'_d.tar ']);
unix(['cd ~/CS/tmp;mv ~/CS/BAC/',outFileDirName,'_t.tar .;',' tar xf ',outFileDirName,'_t.tar ']);
