if 1==2
  clear
  bcs_data_path = '~/CS/BAC/'
  s2 = load([bcs_data_path,'Niv/SR50/sample2_dat']);
  s3 = load([bcs_data_path,'Niv/SR50/sample3_dat']);
   
  for i=1:6
    subplot(2,3,i);
    plot(s2.F{i},s3.F{i},'.')
    xlabel('F in 2');ylabel('F in 3')
    ay = get(gca,'ylim');
    ax = get(gca,'xlim');
    mx= max(ay,ax);
    line([0,mx],[0 mx]);axis('equal')
    title(i)
  end
  
  print -dpdf ~/CS/BAC/Niv/SR50/FbyRegions

  
end


clear
bcs_data_path = '~/CS/BAC/'
s3=  load([bcs_data_path,'Niv/SR50/sample3_res_indices_moreThan1e3_Jan2013']);
s2=  load([bcs_data_path,'Niv/SR50/sample2_res_indices_moreThan1e3_Jan2013']);

%load([bcs_data_path,'full16S/data_primersAmit_0MM_readLength_26'])


clf
plot(s2.res,s3.res,'.');
ay = get(gca,'ylim');
ax = get(gca,'xlim');
mx= max(ay,ax);
line([0,mx],[0 mx])


xlabel('sample 2');ylabel('sample 3');title('bacteria found in two technical repeats. read length 26')
set(gca,'fontweight','bold','fontsize',15)
print -dpdf ~/CS/BAC/Niv/SR50/comparison_sample23




%xlabel('2');ylabel('3');title('regions 2,3,4,5,6');
%print -dpdf ~/CS/BAC/Niv/SR50/comparison_regions23456

%xlabel('2');ylabel('3');title('all regions');
%print -dpdf ~/CS/BAC/Niv/SR50/comparison_allRegions

%xlabel('2');ylabel('3');title('regions 1,2,3,5,6');
%print -dpdf ~/CS/BAC/Niv/SR50/comparison_regions12356


a = find(s2.res>0.02 & s3.res<0.013 | s2.res<0.02 & s3.res>0.013)

%a = find(s2.res>0.03 & s3.res<0.03 | s2.res<0.03 & s3.res>0.03)

for i=1:length(a3)
  clear d
  d = NaN*ones(1,length(s2.found_bact_in_26));
  for j=1:length(s2.found_bact_in_26)
    d(j) = length(find(seq_readLength26(a3(i),:)-seq_readLength26(s2.found_bact_in_26(j),:)));
  end
  dist_a3(i) = min(d);
end

for i=1:length(a2)
  clear d
  d = NaN*ones(1,length(s3.found_bact_in_26));
  for j=1:length(s3.found_bact_in_26)
    d(j) = length(find(seq_readLength26(a2(i),:)-seq_readLength26(s3.found_bact_in_26(j),:)));
  end
  dist_a2(i) = min(d);
end

hist([dist_a2,dist_a3],100)
title('distance between non-overlapping solutions samples 2 and 3')
print -dpdf ~/CS/BAC/Niv/SR50/distanceOfNonOverlapping

204721
