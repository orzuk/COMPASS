function [strinfilename,red]=prepareReadsForRandomSequences(Sequence1,seed_in,Nbac_in_mixture,Nbac_in_mixture)


disp('prepares reads for random sequences - for the the infinite case')
%keyboard

if isempty(seed_in)
  outseed = sum(100*clock);
else
  outseed=seed_in;  
end
rand('seed',outseed);
%keyboard

w=1./([1:Nbac_in_mixture].^npower);%setting the frequencies as the power law
w=round(w/sum(w)*Nreads);%this is how many read  we will have from each bacteria in the mix
if bac_dist_flag==0
  disp('random selection')  
  ind_bac_in_mix=randi([1,N],Nbac_in_mixture,1);
else
  disp('close bacteria')
  tmp=randi([1+Nbac_in_mixture,length(indsimbac_vec)-Nbac_in_mixture],1,1);
  ind_bac_in_mix=indsimbac_vec(tmp-ceil(Nbac_in_mixture/2):tmp+ceil(Nbac_in_mixture/2)-1);
end

ind_bac_in_mix=ind_bac_in_mix(randperm(Nbac_in_mixture));
%keyboard

if exist('ind_bac_in_mix_input')
  ind_bac_in_mix = ind_bac_in_mix_input.nonZeroBAC;
  allBAC  = ind_bac_in_mix_input.allBAC;
  disp('entered ind_bac_in_mix - changed this')
  pause
end


%keyboard
correctWeight = zeros(N,1);
%11,correctWeight(251001:252000) = 0;
correctWeight(ind_bac_in_mix ) = w./sum(w);


npowerName = num2str(npower);
npowerName(find(npowerName=='.')) = [];

basicStrintFileName = [outfilename,'_Nbacmix_',num2str(Nbac_in_mixture),'_Nread_',num2str(Nreads)...
    ,'_Readlen_',num2str(readlen),'_Npower_',npowerName,'_bacdistflag_',num2str(bac_dist_flag)];

%strinfilename=[basicStrintFileName,'_reads.mat'];
strinfilename = [basicStrintFileName,'_NoReads.mat'];




if isempty(seed_in) | saveToFileFlag==1 % if no seed - saves to file
  %save(strinfilename,'Nbac_in_mixture', 'Header_mix','ind_bac_in_mix', 'w', 'outseed','correctWeight','Nreads','readlen','npower','bac_dist_flag');
  save(strinfilename,'Nbac_in_mixture', 'ind_bac_in_mix', 'w', 'outseed','correctWeight','Nreads','readlen','npower','bac_dist_flag');
end




if createReadsFlag==0
  red = [];
  return
end



%keyboard
load(basicSeqKey,'positionInPart')
Sequence_mix = cell(length(ind_bac_in_mix),1);
for i=1:length(ind_bac_in_mix)
  clear seq_*
  load([basicSeqNameDir,'seq_part_',num2str(positionInPart(ind_bac_in_mix(i)))],['seq_',num2str(ind_bac_in_mix(i))]);
  ww = ['Sequence_mix{i} = seq_',num2str(ind_bac_in_mix(i)),'{1};'];
  eval(ww);
  
  clear head_*
  load([basicSeqNameDir,'head_part_',num2str(positionInPart(ind_bac_in_mix(i)))],['head_',num2str(ind_bac_in_mix(i))]);
  ww = ['Header_mix{i} = head_',num2str(ind_bac_in_mix(i)),'{1};'];
  eval(ww);
end
clear seq_* positionInPart head_*

%keyboard
len_mix=cellfun(@length,Sequence_mix);
bac_read_vec=zeros(Nreads,1);
bac_read_startvec=zeros(Nreads,1);


k=1;
for i=1:Nbac_in_mixture
    bac_read_vec(k:k+w(i)-1)=i;
    bac_read_startvec(k:k+w(i)-1)=randi([1,len_mix(i)-readlen+1],w(i),1);
    k=k+w(i);
end
red=repmat('n',Nreads,readlen);
k=1;
for i=1:Nbac_in_mixture
    for j=1:w(i)
        red(k+j-1,:)=Sequence_mix{i}(bac_read_startvec(k+j-1):bac_read_startvec(k+j-1)+readlen-1);
    end
    k=k+w(i);
end

tmp=randperm(Nreads);
red=red(tmp,:);
bac_read_vec=bac_read_vec(tmp);

%keyboard
%correctWeight = zeros(1,455055);
%correctWeight(ind_bac_in_mix) = w./sum(w);
%correctWeight = correctWeight';
%keyboard