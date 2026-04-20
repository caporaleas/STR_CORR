function [corr_acrossSubjs_allGM,corr_acrossSubjs_allWM,SE_R1_GM] = STR_CORR_12_whole_GM(initpath,outpath,parameters,V_R1_maps_filtered,V_diff_maps_filtered,PVE_masks,AtlasGM_masks,AtlasWM_masks,confidence_bounds)

%% Calculates medians of the parameters in the entire GM (masked with PVE and the Atlas)
%calculates median of R1 and diff parameter within each ROI and for each
%subj
cd(initpath)
folders = dir('pil*');
N_subj = 20;  % Number of healthy subjects
N_par = length(parameters); % Number of metrics

corr_acrossSubjs_allGM = zeros(N_par-1,2);
corr_acrossSubjs_allWM = zeros(N_par-1,2);

median_R1_GM = zeros(N_subj,1);
median_R1_WM = zeros(N_subj,1);
median_diff_GM = zeros(N_subj,1,N_par-1);
median_diff_WM = zeros(N_subj,1,N_par-1);
mad_R1_GM = zeros(N_subj,1);
mad_R1_WM = zeros(N_subj,1);
mad_diff_GM = zeros(N_subj,1,N_par-1);
mad_diff_WM = zeros(N_subj,1,N_par-1);
SE_R1_GM = zeros(N_subj,1);
SE_R1_WM = zeros(N_subj,1);
SE_diff_GM = zeros(N_subj,1,N_par-1);
SE_diff_WM = zeros(N_subj,1,N_par-1);


tmp0 = AtlasGM_masks{1,1};
tmp00 = zeros(size(tmp0));

for rr = 1 : length(AtlasGM_masks)
    tmp00 = tmp00 + AtlasGM_masks{rr,1};
end

AtlasGM_entiremask = tmp00;

tmp0 = AtlasWM_masks{1,1};
tmp00 = zeros(size(tmp0));

for rr = 1 : length(AtlasWM_masks)
    tmp00 = tmp00 + AtlasWM_masks{rr,1};
end

AtlasWM_entiremask = tmp00;

% GM and WM ROIs
disp('Calculate medians for the entire GM and WM tissues')
for ii = 1 : N_subj 
    sprintf('Processing subj %d',ii)

    disp('Calculating R1 medians')

    tmp = V_R1_maps_filtered{ii,1}.*PVE_masks{ii,1}.*AtlasGM_entiremask;
    V_R1_masked_arr=tmp(tmp~=0);
    median_R1_GM(ii,1) = median(V_R1_masked_arr);
    mad_R1_GM(ii,1) = mad(median_R1_GM,1,1);
    disp(mad_R1_GM(ii,1))
    SE_R1_GM(ii,1) = mad_R1_GM(ii,1).*1.4826./sqrt(length(V_R1_masked_arr));

    tmp = V_R1_maps_filtered{ii,1}.*PVE_masks{ii,2}.*AtlasWM_entiremask;
    V_R1_masked_arr=tmp(tmp~=0);
    median_R1_WM(ii,1) = median(V_R1_masked_arr);
    mad_R1_WM(ii,1) = mad(median_R1_WM,1,1);
    SE_R1_WM(ii,1) = mad_R1_WM(ii,1).*1.4826./sqrt(length(V_R1_masked_arr));

    disp('Calculating diffusion metrics medians')
    for pp = 1 : N_par-1

        tmp = V_diff_maps_filtered{ii,pp}.*PVE_masks{ii,1}.*AtlasGM_entiremask;
        V_diff_masked_arr=tmp(tmp~=0);
        median_diff_GM(ii,1,pp) = median(V_diff_masked_arr);
        mad_diff_GM(ii,1,pp) = mad(median_diff_GM(ii,1,pp),1,1);
        SE_diff_GM(ii,1,pp) = mad_diff_GM(ii,1,pp).*1.4826./sqrt(length(V_diff_masked_arr));

        tmp = V_diff_maps_filtered{ii,pp}.*PVE_masks{ii,2}.*AtlasWM_entiremask;
        V_diff_masked_arr=tmp(tmp~=0);
        median_diff_WM(ii,1,pp) = median(V_diff_masked_arr);
        mad_diff_WM(ii,1,pp) = mad(median_diff_WM(ii,1,pp),1,1);
        SE_diff_WM(ii,1,pp) = mad_diff_WM(ii,1,pp).*1.4826./sqrt(length(V_diff_masked_arr));

    end
end

correlations = zeros(N_par-1,1);
rawpvalues = zeros(N_par-1,1);


disp('Calculate correlations bw R1 and diff metrics across subjs')
plotname_GM = 'Scatterplot_GM_diffpar_vs_R1_wholeGM_acrossSubjs';
h = figure;
for pp = 1 : N_par-1
    
    subplot(2,4,pp) 

    % median across subjects for each region
       
    [r,p] = corrcoef(median_R1_GM, median_diff_GM(:,1,pp), 'rows','complete');

    correlations(pp,1) = r(2);
    rawpvalues(pp,1) = p(2);

    corr_acrossSubjs_allGM(pp,1) = r(2);
    corr_acrossSubjs_allGM(pp,2) = p(2);

    corr_coef = round(r(2),3,'significant');
    p_value = round(p(2),3,'significant');

    corr_coef_str_GM = num2str(corr_coef);
    p_value_str_GM = num2str(p_value);
    
    h.Color = [1 1 1];
    h.Position = [1 1 1200 800];

    s = errorbar(median_R1_GM, median_diff_GM(:,1,pp), SE_diff_GM(:,1,pp), SE_diff_GM(:,1,pp), SE_R1_GM, SE_R1_GM,'o');
    s.LineWidth = 0.6;
    s.MarkerEdgeColor = 'b';
    s.MarkerFaceColor = [0.8 0.8 0.8];
    s.Color = 'b';
    hold on
    plot(median_R1_GM, median_diff_GM(:,1,pp),'bo');

    if confidence_bounds == 1
        % confidence bounds
        x = median_R1_GM;
        y = median_diff_GM(:,1,pp);
        [p,S] = polyfit(x,y,1);
        % Ensure the confidence interval covers the full x-axis
        xLimits = xlim;  % Get current x-axis limits
        xv = linspace(xLimits(1), xLimits(2), 150);  % Generate x values over full x-axis
        [y_ext,delta] = polyconf(p,xv,S);
        hold on
        %plot(xv, y_ext, '-b')
        patch([xv fliplr(xv)], [(y_ext+delta) fliplr((y_ext-delta))], 'b', 'FaceAlpha',0.1, 'EdgeColor','none')
    end

    hold on

    kk = lsline;
    kk.Color = 'b'; kk.LineWidth = 2;
    ylabel(parameters{1,pp+1},'FontSize',14,'FontWeight','Bold');
    xlabel("R1 (s^{-1})",'FontSize',14,'FontWeight','Bold');
    %xlim([0.5 1]); ylim([0 0.5]);
 
    title({strcat('r_{GM}=',corr_coef_str_GM,'; p=',p_value_str_GM)});

    ax = gca;
    ax.Box = 'off'; % Remove the box around the plot
    ax.FontSize = 12;

end
saveas(gcf,fullfile(outpath,plotname_GM),'fig');

% FDR correction using Benjamini-Hochberg method
if exist('mafdr', 'file') % Check if mafdr is available
    adjpvalues = mafdr(rawpvalues, 'BHFDR', true);
else
    % Manual FDR correction if mafdr is not available
    [sortedP, sortIdx] = sort(rawpvalues);
    numTests = length(pValues);
    adjustedP = sortedP .* numTests ./ (1:numTests)';
    adjustedP(adjustedP > 1) = 1; % Ensure p-values are not >1
    
    % Reorder back to original indices
    adjpvalues = zeros(size(rawpvalues));
    adjpvalues(sortIdx) = adjustedP;
end

length(correlations)
resultsTable = table(parameters(2:end)',correlations, rawpvalues, adjpvalues, ...
                     'VariableNames', {'Parameter','Whole_GM_Correlation', 'Raw_pValue', 'Adjusted_pValue'});
disp(resultsTable);

end