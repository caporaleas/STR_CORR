function [V_pve_masks] = STR_CORR_02_getting_pve(initpath,outpath,PVEerosion)
%STR_CORR_02_getting_pve gets partial volume estimate maps and binarize them based on the absolute maximum
% among the three tissues (GM, WM, CSF)
% Directory with subject data
cd(initpath)
folders = dir('pil*');
N_subj = length(folders);  % Number of subjects

% cell initialization
V_pve_maps = cell(N_subj,3);
V_pve_masks = cell(N_subj,3);

% getting WM, GM, CSF pve from T1 segmentation (to MNI)
for ii = 1 : N_subj  
    disp('pve')
    disp(ii)

    img_path_GMpve = strcat('/certainpath/Charmed/Dati_GM_mask/',folders(ii).name,'/',folders(ii).name,'_pve_1_MNI_new.nii.gz');
    nii = load_untouch_nii(img_path_GMpve);
    V_GM_tot = double(nii.img);
    V_pve_maps{ii,1} = V_GM_tot;
    
    img_path_WMpve = strcat('/certainpath/Charmed/Dati_GM_mask/',folders(ii).name,'/',folders(ii).name,'_pve_2_MNI_new.nii.gz');
    nii = load_untouch_nii(img_path_WMpve);
    V_WM_tot = double(nii.img);
    V_pve_maps{ii,2} = V_WM_tot;
    
    img_path_CSFpve = strcat('/certainpath/Charmed/Dati_GM_mask/',folders(ii).name,'/',folders(ii).name,'_pve_0_MNI_new.nii.gz');
    nii = load_untouch_nii(img_path_CSFpve);
    V_CSF_tot = double(nii.img);
    V_pve_maps{ii,3} = V_CSF_tot;
end

for ii = 1 : N_subj

    disp('GM, WM and CSF pve')
    disp(ii)

    GMpve = V_pve_maps{ii,1};
    WMpve = V_pve_maps{ii,2};
    CSFpve = V_pve_maps{ii,3};

    GMpve_mask = GMpve;
    GMpve_mask(GMpve_mask>0)=1;
    WMpve_mask = WMpve;
    WMpve_mask(WMpve_mask>0)=1;
    CSFpve_mask = CSFpve;
    CSFpve_mask(CSFpve_mask>0)=1;

    for zz = 1 : size(GMpve,3)
        for yy = 1 : size(GMpve,2)
            for xx = 1 : size(GMpve,1)
                pve_arr(1,1)=GMpve(xx,yy,zz);
                pve_arr(1,2)=WMpve(xx,yy,zz);
                pve_arr(1,3)=CSFpve(xx,yy,zz);
                [m,n] = find(pve_arr==max(pve_arr));
                if length(m)==1
                    if n == 1
                        WMpve_mask(xx,yy,zz)=0;
                        CSFpve_mask(xx,yy,zz)=0;
                    elseif n == 2
                        GMpve_mask(xx,yy,zz)=0;
                        CSFpve_mask(xx,yy,zz)=0;
                    elseif n == 3
                        WMpve_mask(xx,yy,zz)=0;
                        GMpve_mask(xx,yy,zz)=0;
                    end
                end
            end
        end
    end
    if PVEerosion == 1
        for ss = 1 : size(GMpve_mask,3)
            tmp = GMpve_mask(:,:,ss);
            bw1 = edge(tmp);
            tmp = tmp-bw1;
            tmp(tmp<0) = 0;
            GMpve_mask(:,:,ss) = tmp;
        end
        V_pve_masks{ii,1} = GMpve_mask;

        for ss = 1 : size(WMpve_mask,3)
            tmp = WMpve_mask(:,:,ss);
            bw2 = edge(tmp);
            tmp = tmp-bw2;
            tmp(tmp<0) = 0;
            WMpve_mask(:,:,ss) = tmp;
        end
        V_pve_masks{ii,2} = WMpve_mask;

        for ss = 1 : size(CSFpve_mask,3)
            tmp = CSFpve_mask(:,:,ss);
            bw3 = edge(tmp);
            tmp = tmp-bw3;
            tmp(tmp<0) = 0;
            CSFpve_mask(:,:,ss) = tmp;
        end
        V_pve_masks{ii,3} = CSFpve_mask;
        
        disp('Erosion of pve masks applied')
        figname = 'eroded_pve_maps';
    else
        V_pve_masks{ii,1} = GMpve_mask;
        V_pve_masks{ii,2} = WMpve_mask;
        V_pve_masks{ii,3} = CSFpve_mask;
        disp('No erosion applied on PVE masks')
        figname = 'pve_maps';
    end
end

% PVE map showing WM, GM and CSF
N_slice = [32 38 48 58 68];

h = figure;
for ss = 1 : length(N_slice)
    slice = N_slice(ss);
    subplot(1,length(N_slice),ss)
    tmp1 = V_pve_masks{1,1}(:,:,slice);
    tmp2 = V_pve_masks{1,2}(:,:,slice);
    tmp3 = V_pve_masks{1,3}(:,:,slice);
    tmp = tmp1 + tmp2.*2 + tmp3.*3;
    imshow(rot90(tmp),[]); colormap('jet'); c = colorbar; c.Location = 'southoutside'; c.Color = [1 1 1];
end
h.Color = [0 0 0.513725];
saveas(h,fullfile(outpath,figname),'fig');

close all

end
