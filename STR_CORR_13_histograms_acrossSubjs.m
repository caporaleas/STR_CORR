function [corr_for_each_ROI_GM,corr_for_each_ROI_WM] = STR_CORR_13_histograms_acrossSubjs(median_R1_GM,median_R1_WM,median_diff_GM,median_diff_WM,parameters,outpath)
% calculates correlation coefficient matrix computing the median of the
% medians across ROIs and creates histograms of r

N_par = length(parameters);
N_GM_ROIs = size(median_R1_GM,2);
N_WM_ROIs = size(median_R1_WM,2);
N_subjs = size(median_R1_GM,1);


plotname_GM = 'Histograms_GM_diffpar_vs_R1_acrossSubjs';
plotname_WM = 'Histograms_WM_diffpar_vs_R1_acrossSubjs';

corr_for_each_ROI_GM = zeros(N_subjs,N_par-1,2);
corr_for_each_ROI_WM = zeros(N_subjs,N_par-1,2);

for ii = 1 : N_GM_ROIs
    for pp = 1 : N_par-1
        
        [r,p] = corrcoef(median_R1_GM(:,ii), median_diff_GM(:,ii,pp), 'rows','complete');
        
        corr_for_each_ROI_GM(ii,pp,1) = r(2);
        corr_for_each_ROI_GM(ii,pp,2) = p(2);
    end
end
for ii = 1 : N_WM_ROIs
    for pp = 1 : N_par-1

        [r,p] = corrcoef(median_R1_WM(:,ii), median_diff_WM(:,ii,pp), 'rows','complete');

        corr_for_each_ROI_WM(ii,pp,1) = r(2);
        corr_for_each_ROI_WM(ii,pp,2) = p(2);

    end
end

h = figure;
h.Color = [1 1 1];
for pp = 1 : N_par-1
    subplot(2,4,pp)
    hh = histogram(corr_for_each_ROI_GM(:,pp,1),10); title(parameters{1,pp+1});
    hh.FaceColor = 'b';
    hh.EdgeColor = 'k';
end

saveas(gcf,fullfile(outpath,plotname_GM),'fig');

h = figure;
h.Color = [1 1 1];
for pp = 1 : N_par-1
    subplot(2,4,pp)
    hh = histogram(corr_for_each_ROI_WM(:,pp,1),10); title(parameters{1,pp+1});
    hh.FaceColor = 'r';
    hh.EdgeColor = 'k';
end

saveas(gcf,fullfile(outpath,plotname_WM),'fig');



end
