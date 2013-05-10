function [PrimerReadBySpeciesMat Primer_kmers_packed bactInRegion]=oneRegionDataForConcatenatedData(reg,readLength,SS)

%keyboard

bactInRegion = find(SS.amplifiedRegion(:,reg));
basicOffSet = readLength*6;
p_position_inds_f = (reg-1)*readLength+1;
p_position_inds_r = basicOffSet+(reg-1)*readLength+1;

p_position_inds = [p_position_inds_f*ones(1,length(bactInRegion)),p_position_inds_r*ones(1,length(bactInRegion))];

bact_written_twice = [bactInRegion;bactInRegion];

[PrimerReadBySpeciesMat Primer_kmers_packed] = BuildMixingMatrixFromSequencesRegions( ...
    readLength, SS.Sequence_packed, SS.len_uni, 64, bact_written_twice,p_position_inds); % ma

if ~isempty(intersect(int2nt(unpack_seqs(Primer_kmers_packed,readLength,64)),char('A')*ones(1,readLength),'rows'))
  disp('problem oneRegionDataForConcatenatedData.m');
  pause
end
