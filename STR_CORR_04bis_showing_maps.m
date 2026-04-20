function done3 = STR_CORR_04bis_showing_maps(inpath,outpath)
%PLOTS PARAMETRIC MAPS IN MNI (AVERAGE ACROSS SUBJECTS)

cd(inpath)
cd('Dati_R1')

nii = load_untouch_nii('R12std_SubjMean_ero.nii.gz');
R12std_ero = double(nii.img);

cd(inpath)
cd('Dati_SANDI')
nii = load_untouch_nii('Rsoma2std_SubjMean_ero.nii.gz');
Rsoma2std_ero = double(nii.img);
nii = load_untouch_nii('fsoma2std_SubjMean_ero.nii.gz');
fsoma2std_ero = double(nii.img);
nii = load_untouch_nii('fneurite2std_SubjMean_ero.nii.gz');
fneurite2std_ero = double(nii.img);
nii = load_untouch_nii('fextra2std_SubjMean_ero.nii.gz');
fextra2std_ero = double(nii.img);
nii = load_untouch_nii('De2std_SubjMean_ero.nii.gz');
De2std_ero = double(nii.img);
nii = load_untouch_nii('Din2std_SubjMean_ero.nii.gz');
Din2std_ero = double(nii.img);

N_slice = [32 38 48 58 68];
%disp(class(V_R1_maps));  % Should display 'cell'
%disp(class(V_diff_maps)); % Should display 'cell'

h = figure;
for ss = 1 : length(N_slice)
    slice = N_slice(ss);
    subplot(1,length(N_slice),ss)
    tmp = R12std_ero(:,:,slice);
    imshow(rot90(tmp),[0 1.2]); colormap('copper'); c = colorbar; c.Location = 'southoutside'; c.Color = [1 1 1];
end
h.Color = [0 0 0];
saveas(h,fullfile(outpath,'R1_maps_MNI'),'fig');

h = figure;
for ss = 1 : length(N_slice)
    slice = N_slice(ss);
    subplot(1,length(N_slice),ss)
    tmp = fsoma2std_ero(:,:,slice);
    imshow(rot90(tmp),[0 0.6]); colormap('hot'); c = colorbar; c.Location = 'southoutside'; c.Color = [1 1 1];
end
h.Color = [0 0 0];
saveas(h,fullfile(outpath,'fsoma_maps_MNI'),'fig');

h = figure;
for ss = 1 : length(N_slice)
    slice = N_slice(ss);
    subplot(1,length(N_slice),ss)
    tmp = fneurite2std_ero(:,:,slice);
    imshow(rot90(tmp),[0 0.6]); colormap('hot'); c = colorbar; c.Location = 'southoutside'; c.Color = [1 1 1];
end
h.Color = [0 0 0];
saveas(h,fullfile(outpath,'fneurite_maps_MNI'),'fig');

h = figure;
for ss = 1 : length(N_slice)
    slice = N_slice(ss);
    subplot(1,length(N_slice),ss)
    tmp = fextra2std_ero(:,:,slice);
    imshow(rot90(tmp),[0 0.8]); colormap('hot'); c = colorbar; c.Location = 'southoutside'; c.Color = [1 1 1];
end
h.Color = [0 0 0];
saveas(h,fullfile(outpath,'fextra_maps_MNI'),'fig');

h = figure;
for ss = 1 : length(N_slice)
    slice = N_slice(ss);
    subplot(1,length(N_slice),ss)
    tmp = De2std_ero(:,:,slice);
    imshow(rot90(tmp),[0 3]); colormap('hot'); c = colorbar; c.Location = 'southoutside'; c.Color = [1 1 1];
end
h.Color = [0 0 0];
saveas(h,fullfile(outpath,'De_maps_MNI'),'fig');

h = figure;
for ss = 1 : length(N_slice)
    slice = N_slice(ss);
    subplot(1,length(N_slice),ss)
    tmp = Din2std_ero(:,:,slice);
    imshow(rot90(tmp),[1.5 2.5]); colormap('hot'); c = colorbar; c.Location = 'southoutside'; c.Color = [1 1 1];
end
h.Color = [0 0 0];
saveas(h,fullfile(outpath,'Din_maps_MNI'),'fig');

h = figure;
for ss = 1 : length(N_slice)
    slice = N_slice(ss);
    subplot(1,length(N_slice),ss)
    tmp = Rsoma2std_ero(:,:,slice);
    imshow(rot90(tmp),[10 15]); colormap('hot'); c = colorbar; c.Location = 'southoutside'; c.Color = [1 1 1];
end
h.Color = [0 0 0];
saveas(h,fullfile(outpath,'Rsoma_maps_MNI'),'fig');

close all




done3 = 'Parametric maps in MNI saved into figure';
disp(done3);
end