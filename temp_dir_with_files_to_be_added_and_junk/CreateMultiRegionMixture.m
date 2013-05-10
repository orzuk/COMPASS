% numbact - number of bacteria in the mixture (ignored if a list of bactid is given)
% bactid - list of IDs of the bacteria or empty to choose random
% bactfreq - frequencies of the bacteria or empty to choose random
% distrib (ignored if bactfreq is non-empty) - two options -
% - 0 constant frequency 1/numbact
% - 1 powerlaw -1 frequency (not implemented yet)
% numreadsperregion - number of reads to generate for each region
% readlen - for each read
% regions - a list of amplified regions (start,end for each bacteria for
% each region) - an array (of regions) containing pos which is N*2 list
% (sp from Data/SequenceDB/primers_position_few_regions.forsimulation.nonambig.mat)
% seqdb - the bacteria sequence database from which to generate the reads
% (use the non-redundant database 410849 seq db)
%
% returns
% reads - the reads matrix (numreads x 4 uint64)
% bact,freq - the bacteria in the mixture and their frequency (indices from seqdb)
function [reads,bact,freq]=CreateMultiRegionMixture(numbact,bactid,bactfreq,distrib,numreadsperregion,readlen,regions,seqdb)

numseqdb=size(seqdb,2);

% generate the bacteria id
if (isempty(bactid))
    rbperm=randperm(numseqdb);
    bact=rbperm(1:numbact);
else
    bact=bactid;
    numbact=length(bact);
end

% generate bacteria frequencies
if (isempty(bactfreq))
    if (distrib==0)
        freq=ones(numbact,1);
    else
        freq=rand(numbact,1);
        disp('option non supported');
    end
    freq=freq/sum(freq);
else
    freq=bactfreq;
end

numregions=length(regions);

reads=zeros(numreadsperregion*numregions,ceil(readlen/32),'uint64');
%fromb=zeros(numreadsperregion*numregions,1);
%fromp=zeros(numreadsperregion*numregions,1);
seqs=seqdb(bact);

oread=1;
for cregion=1:numregions
    % create list of all amplified bacteria for this region
    cregdat=regions(cregion).pos;
    bactpresent=((cregdat(bact,1)>0)&(cregdat(bact,2)>0));
    abact=bact(bactpresent);
    numabact=length(abact);
    cfreq=cumsum(freq(bactpresent))./sum(freq(bactpresent));
    if (~isempty(cfreq))
        for cread=1:numreadsperregion
            % select bacteria to get sequence from
            cnum=rand(1,1);
            cbactid=find(cnum<=cfreq,1,'first');
            cbact=abact(cbactid);
            % select read direction
            cdir=randi(2,1);
            if (cdir==1)
                readstart=cregdat(cbact,1);
            else
                readstart=cregdat(cbact,2)+1-readlen;
            end
            
            cseq=seqdb{cbact};
            reads(oread,:)=pack_seqs(cseq(readstart:readstart+readlen-1),64);
 %           fromb(oread)=cbact;
 %           fromp(oread)=readstart;
            oread=oread+1;
        end
    else
        disp('no bacteria amplified');
    end
end
