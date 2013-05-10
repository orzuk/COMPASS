load('~/CS/BAC/primers_startAs750_length350/sample_read_qual') 
I = [11 3 6 10 8 1 4 7 9 2 12 5]

fid = fopen('~/CS/BAC/12Samples/data454.txt','w');
for i=I
  clear  uniqueReads uniqueReads_length
  load(['~/CS/BAC/12Samples/454/data/reads_samples_',sample_read_qual{i,1},'_andMatch_basedFullDataBaseOf350WithAmbiguous'])
  %fprintf(fid,'%s\t%f1.0\n',sample_read_qual{i,1},sum(uniqueReads_length));
  fprintf(fid,'%s\t ~%1.0e\n',sample_read_qual{i,1},sum(uniqueReads_length));
end
fclose(fid);

clear
samples = {'S1','S2','S3','S4','M1','M2','M3','M4','O7','O10','S7','S10'}

fid = fopen('~/CS/BAC/12Samples/dataSolexa.txt','w');
for i=1:length(samples)
  load(['~/CS/BAC/12Samples/Solexa/data/primer_tail/illumina_reads100_',samples{i},'_uni_primer_tail_noR']);
  w = ['a = sum(freq_uni_reads100_',samples{i},',1);'];
  eval(w);
  
  fprintf(fid,'%s\t ~%1.0e\n',samples{i},round(a/10^6)*10^6);
end
fclose(fid);


