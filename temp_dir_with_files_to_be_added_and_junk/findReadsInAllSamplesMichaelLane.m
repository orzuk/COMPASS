clear
bcs_data_path = '~/CS/BAC/'

for sampleNum=1:20
  for i=1:6
    tmp = load([bcs_data_path,'MichaelLane/data/sample_',num2str(sampleNum),'_region_',num2str(i),'_unireads.mat']);
    reads{sampleNum}{i}.uniqueReads = pack_seqs(tmp.readsuni,64);
    reads{sampleNum}{i}.uniqueReads_length = tmp.frequni;
  end
end

clear keep
for i=1:6
  keep{i} = reads{1}{i}.uniqueReads;
end
  
for i=1:20
  for j=i+1:20
    for k=1:6
      [junk,i1,i2] = intersect(reads{i}{k}.uniqueReads,reads{j}{k}.uniqueReads,'rows');
      keep{k} = intersect(keep{k},junk,'rows');
    end
  end
end

% how many of them in each sample

freq = cell(6,1);
for i=1:6
  freq{i} = zeros(size(keep{i},1),20);
end
for sampleNum=1:20
  for k=1:6
    [junk,i1,i2] = intersect(reads{sampleNum}{k}.uniqueReads,keep{k},'rows');
    if length(i2)~=size(keep{k},1)
      disp('problem')
      pause
    end
    
    freq{k}(i2,sampleNum) = reads{sampleNum}{k}.uniqueReads_length(i1)/sum(reads{sampleNum}{k}.uniqueReads_length);
    
  end
end



for i=1:6
  figure(i)
  for j=1:20
    subplot(4,5,j);
    plot(freq{i}(:,j));
    title([i,j])
  end
end

g1 = [1:7 9:15];
g2 = [16:20];
for i=1:6
  f1 = freq{i}(:,g1);
  f2 = freq{i}(:,g2);
  [x1,y1] = find(f1>0.05);
  [x2,y2] = find(f2>0.05);
  length(intersect(x1,x2))
end

I looked at the reads that appear in all samples (in each region).
although there is a large amount of them - it never happens that a high frequency read in the Michael samples appear as high in Peter's samples - switches happen in 0.01%.
Namely a read may appear as `% in the other group.


clear totalOverlap
for i=1:6
  totalOverlap(i,:) = sum(freq{i},1);
end

% are they mapped to something else?



figure(1);clf
for i=1:6
  subplot(2,3,i)
  imagesc(freq{i},[0 0.4])
  %colorbar
  title(i)
end
figure(1)
print ~/CS/BAC/MichaelLane/results/overlappingReads

figure(2)
imagesc(totalOverlap)
colorbar
title('total fraction of overlapping reads')
print('-dpdf','~/CS/BAC/MichaelLane/results/totalOverlapAmongRegions')  


% what are the bacteria which correspond to these high frequency?
for sampleNum=1:20
  for k=1:6
    [junk,i1,i2] = intersect(reads{sampleNum}{k}.uniqueReads,keep{k},'rows');
    if length(i2)~=size(keep{k},1)
      disp('problem')
      pause
    end
    
    freq{k}(i2,sampleNum) = reads{sampleNum}{k}.uniqueReads_length(i1)/sum(reads{sampleNum}{k}.uniqueReads_length);
    
  end
end
