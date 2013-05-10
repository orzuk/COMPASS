function [PrimerReadBySpeciesMat Primer_kmers_packed bactInRegion]=oneRegionDataForPE(reg,readLength,SS)
disp('is the readLength multiplied by 2?')
%keyboard

bactInRegion = find(SS.amplifiedRegion(:,reg));
p_position_inds_f = (reg-1)*readLength+1;

p_position_inds = [p_position_inds_f*ones(1,length(bactInRegion))];

[PrimerReadBySpeciesMat Primer_kmers_packed] = BuildMixingMatrixFromSequencesRegions( ...
    readLength, SS.Sequence_packed, SS.len_uni, [], bactInRegion,p_position_inds); % ma

if ~isempty(intersect(int2nt(unpack_seqs(Primer_kmers_packed,readLength,64)),char('A')*ones(1,readLength),'rows'))
  disp('problem oneRegionDataForPE.m');
  pause
end
