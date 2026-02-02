function Barrett_analysis_no_correction_cluster(ts,sesidx,savedir,datastream,iter)

%  kept exact CopBET to instead take single TS and opperate within parfor loop,
% saving outputs including disp along the way.

if isempty(iter)
    iter=1000;
end

diary([savedir,'/outlogs/outlog_',datastream,'_ses',num2str(sesidx),'.txt']);


Ct2 = DCC(ts-mean(ts),'iter',iter);

num_rois = size(Ct2,1);
variance_out = var(Ct2,[],3);
entropy_out = zeros(num_rois, num_rois);

for i = 1:num_rois
    for ii = i:num_rois
        edge_hist = histcounts(Ct2(i,ii,:),'Normalization','probability');
        entropy_out(i,ii) = nansum(-edge_hist.*log(edge_hist));
    end
end
entropy_out = entropy_out + entropy_out';


%%%%%%%% Checks
% sensible_data_check(variance_out);
entropy_check = entropy_out;
entropy_check(entropy_check==diag(diag(entropy_check)))=1;
sensible_data_check(entropy_check,'entropy matrix');

% check Ct2 is symmetric:
for i = 1:size(Ct2,3)
if norm(Ct2(:,:,i)-Ct2(:,:,i)')>1e-10
    warning('Ct2 matrix not symmetric')
end
end

% check if Ct2 has values outside of correlation coefficient range
if any(Ct2(:)<-1)||any(Ct2(:)>1)
    warning('Not correlation coefficient range')
end

save([savedir,'DCC_',datastream,'_ses',num2str(sesidx),'.mat'],'variance_out','entropy_out')
diary off
end
