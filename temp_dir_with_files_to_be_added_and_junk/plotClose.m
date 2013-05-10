function plotClose(fileName,figNum)

load ~/CS/BAC/dataForSim/sim500CloseDifferentCases_Nbacmix_1000_Nread_2000000_Readlen_50_Npower_05_bacdistflag_1_NoReads
load(fileName)
%figure(figNum)
for i=1:length(found)
  if ~isempty(found{i})
    %hist(abs(found{i}(ind_bac_in_mix)-correctWeight(ind_bac_in_mix)'),100)
    plot(found{i}(ind_bac_in_mix),correctWeight(ind_bac_in_mix)','.');
    axis('equal')
    title(fileName)
  end
end
length(find(abs(found{i}(ind_bac_in_mix)-correctWeight(ind_bac_in_mix)')<0.002))