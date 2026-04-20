%% 3. CALCULATING Pearson's Corr coefficients

L=143; % N_GM_ROIs
img_path_GMatlas='/storage/shared/Atlas/HCPMMP1_on_MNI152_ICBM2009a_nlin_res_on_MNI_T1_2mm.nii.gz';
img_path_WMatlas='/usr/local/fsl_v6062/data/atlases/JHU/JHU-ICBM-labels-2mm.nii.gz';

atlas_GM = load_untouch_nii(img_path_GMatlas);
atlas_GM = double(atlas_GM.img);

r_map_GEASL = atlas; p_map_GEASL = atlas;
r_map_SEASL = atlas; p_map_SEASL = atlas;
r_map_GESE = atlas; p_map_GESE = atlas;
% L = length(atlas_labels_id); 


rPearson_GEASL = zeros(L,1); pPearson_GEASL = zeros(L,1);
rPearson_SEASL = zeros(L,1); pPearson_SEASL = zeros(L,1);
rPearson_GESE = zeros(L,1); pPearson_GESE = zeros(L,1);

for rr = 1 : N_GM_ROIs
    xarr = median_ASLCVR_GM(:,rr);
    yarr = median_GECVR_GM(:,rr);
    [r,p] = corr(xarr,yarr,'rows','complete','Type','Spearman');
    r_map_GEASL(r_map_GEASL==rr) = r;
    p_map_GEASL(p_map_GEASL==rr) = p;
    rPearson_GEASL(rr,1) = r;
    pPearson_GEASL(rr,1) = p;
    
    xarr = median_ASLCVR_GM(:,rr);
    yarr = median_SECVR_GM(:,rr);
    [r,p] = corr(xarr,yarr,'rows','complete','Type','Spearman');
    r_map_SEASL(r_map_SEASL==rr) = r;
    p_map_SEASL(p_map_SEASL==rr) = p;
    rPearson_SEASL(rr,1) = r;
    pPearson_SEASL(rr,1) = p;

    xarr = median_GECVR_GM(:,rr);
    yarr = median_SECVR_GM(:,rr);
    [r,p] = corr(xarr,yarr,'rows','complete','Type','Spearman');
    r_map_GESE(r_map_GESE==rr) = r;
    p_map_GESE(p_map_GESE==rr) = p;
    rPearson_GESE(rr,1) = r;
    pPearson_GESE(rr,1) = p;
end

for rr = 1 : N_WM_ROIs
    xarr = median_ASLCVR_WM(:,rr);
    yarr = median_GECVR_WM(:,rr);
    [r,p] = corr(xarr,yarr,'rows','complete');
    r_map_GEASL(r_map_GEASL==rr+170) = r;
    p_map_GEASL(p_map_GEASL==rr+170) = p;
    rPearson_GEASL(rr+170,1) = r;
    pPearson_GEASL(rr+170,1) = p;
    
    xarr = median_ASLCVR_WM(:,rr);
    yarr = median_SECVR_WM(:,rr);
    [r,p] = corr(xarr,yarr,'rows','complete');
    r_map_SEASL(r_map_SEASL==rr+170) = r;
    p_map_SEASL(p_map_SEASL==rr+170) = p;
    rPearson_SEASL(rr+170,1) = r;
    pPearson_SEASL(rr+170,1) = p;

    xarr = median_GECVR_WM(:,rr);
    yarr = median_SECVR_WM(:,rr);
    [r,p] = corr(xarr,yarr,'rows','complete');
    r_map_GESE(r_map_GESE==rr+170) = r;
    p_map_GESE(p_map_GESE==rr+170) = p;
    rPearson_GESE(rr+170,1) = r;
    pPearson_GESE(rr+170,1) = p;
end