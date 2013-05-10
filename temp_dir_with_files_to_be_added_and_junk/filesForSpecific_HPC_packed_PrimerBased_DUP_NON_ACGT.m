% duplicate the Headers
clear
load ~/CS/BAC/listOfDoubles_primers750
load /u/01/shental/bacteria/bac16s_primers750


%Sequence_packed64(d) = [];
%len_uni(d) = [];
%Header_uni_amp_uni(d) = [];

% must duplicate in the region we amplify in 454 and 

fp_454 = 'CTCCTACGGGAGGCAGCAG';
rp_454 = 'TTGTGCGGGCCCCCGTCAATT'
rc_rp_454 = seqrcomplement(rp_454);

L = zeros(1,length(d));
L454 = L;
for i=1:length(d)
  currSeq = Sequence_uni_amp_uni{d(i)};
  
  j = find(currSeq~='A' & currSeq~='C'  & currSeq~='G'  & currSeq~='T' );
  L(i) = length(j);
  % create the list of of all sequences
  
  %f454 = findstr(currSeq,fp_454);
  %r454 = findstr(currSeq,rc_rp_454);
 
  %if ~isempty(f454) & ~isempty(r454)
  %  currSeq_amp_by_454 = currSeq(f454:r454);
  %  j454 = find(currSeq_amp_by_454~='A' & currSeq_amp_by_454~='C'  & currSeq_amp_by_454~='G'  & currSeq_amp_by_454~='T' );
  
  %  if length(j454)<10
    
  %end
    
  
  %if length(j)<10
  %  basic=listAmbiguousSeq(currSeq);
  
    
  
  %  if ~isempty(f454) & ~isempty(r454)
  %    currSeq_amp_by_454 = currSeq(f454:r454);
  %    basic454 = listAmbiguousSeq(currSeq_amp_by_454);
  %    [size(basic),size(currSeq),size(basic454)]
  %  else
  %    [size(basic),size(currSeq)]
  %  end
  %  pause
  end
end

% if the bacteria is amplified by the 454 also, make sure that the region amplified by 454 is duplicated.

% the first 350 is multiplied



