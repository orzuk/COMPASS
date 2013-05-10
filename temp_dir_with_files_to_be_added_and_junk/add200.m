function add200(fileName)
% add200('/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/VariableReadLength/VariableReadLength_bac_dist_flag_0_Nbac_in_mixture_1000_npower_10_readlen_200_numIter_9_Nreads_1000000_Noise_1_Correction_1')
% for i=[10],add200(['/seq/orzuk2/compressed_sensing/metagenomics/next_gen/CS/BAC/VariableReadLength/VariableReadLength_bac_dist_flag_0_Nbac_in_mixture_1000_npower_10_readlen_200_numIter_',num2str(i),'_Nreads_1000000_Noise_1_Correction_1']);end

load([fileName])
load([fileName,'_L1data'])
load(auxData.reads_data_fileName,'correctWeight','ind_bac_in_mix')
x1_2_nor = addL1toL2_res(uniqueReads,uniqueReads_length,tmpInd,fileNameL1Res,auxData);
CreateMothurDistNoamL1L2(auxData.saveName,tmpInd',x1_2_nor,resCell{end}(:,1),resCell{end}(:,2),correctWeight,auxData);

README = 'save_store_kp{end} is the basis of both results: L2 takes resCell{end}(:,2)> 10^-5 and L1 calculates L1 over save_store_kp{end}'

