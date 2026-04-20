function done = STR_CORR_04_showing_maps(V_R1_maps,V_diff_maps,outpath)
%STR_CORR_04_showing_maps shows the parametric maps 
N_slice = [32 38 48 58 68];
disp(class(V_R1_maps));  % Should display 'cell'
disp(class(V_diff_maps)); % Should display 'cell'

h = figure;
for ss = 1 : length(N_slice)
    slice = N_slice(ss);
    subplot(1,length(N_slice),ss)
    tmp = V_R1_maps{1,1}(:,:,slice);
    tmp = imgaussfilt(tmp,1);
    imshow(rot90(tmp),[0 1.2]); colormap('copper'); c = colorbar; c.Location = 'southoutside'; c.Color = [1 1 1];
end
h.Color = [0 0 0];
saveas(h,fullfile(outpath,'R1_maps'),'fig');

h = figure;
for ss = 1 : length(N_slice)
    slice = N_slice(ss);
    subplot(1,length(N_slice),ss)
    tmp = V_diff_maps{1,1}(:,:,slice);
    imshow(rot90(tmp),[]); colormap('hot'); c = colorbar; c.Location = 'southoutside'; c.Color = [1 1 1];
end
h.Color = [0 0 0];
saveas(h,fullfile(outpath,'FR_maps'),'fig');

h = figure;
for ss = 1 : length(N_slice)
    slice = N_slice(ss);
    subplot(1,length(N_slice),ss)
    tmp = V_diff_maps{1,2}(:,:,slice);
    imshow(rot90(tmp),[]); colormap('hot'); c = colorbar; c.Location = 'southoutside'; c.Color = [1 1 1];
end
h.Color = [0 0 0];
saveas(h,fullfile(outpath,'fneurite_maps'),'fig');

h = figure;
for ss = 1 : length(N_slice)
    slice = N_slice(ss);
    subplot(1,length(N_slice),ss)
    tmp = V_diff_maps{1,3}(:,:,slice);
    imshow(rot90(tmp),[]); colormap('hot'); c = colorbar; c.Location = 'southoutside'; c.Color = [1 1 1];
end
h.Color = [0 0 0];
saveas(h,fullfile(outpath,'fsoma_maps'),'fig');

h = figure;
for ss = 1 : length(N_slice)
    slice = N_slice(ss);
    subplot(1,length(N_slice),ss)
    tmp = V_diff_maps{1,4}(:,:,slice);
    imshow(rot90(tmp),[]); colormap('hot'); c = colorbar; c.Location = 'southoutside'; c.Color = [1 1 1];
end
h.Color = [0 0 0];
saveas(h,fullfile(outpath,'fextra_maps'),'fig');

h = figure;
for ss = 1 : length(N_slice)
    slice = N_slice(ss);
    subplot(1,length(N_slice),ss)
    tmp = V_diff_maps{1,5}(:,:,slice);
    imshow(rot90(tmp),[10 15]); colormap('hot'); c = colorbar; c.Location = 'southoutside'; c.Color = [1 1 1];
end
h.Color = [0 0 0];
saveas(h,fullfile(outpath,'Rsoma_maps'),'fig');

h = figure;
for ss = 1 : length(N_slice)
    slice = N_slice(ss);
    subplot(1,length(N_slice),ss)
    tmp = V_diff_maps{1,6}(:,:,slice);
    imshow(rot90(tmp),[1 2.5]); colormap('hot'); c = colorbar; c.Location = 'southoutside'; c.Color = [1 1 1];
end
h.Color = [0 0 0];
saveas(h,fullfile(outpath,'Din_maps'),'fig');

h = figure;
for ss = 1 : length(N_slice)
    slice = N_slice(ss);
    subplot(1,length(N_slice),ss)
    tmp = V_diff_maps{1,7}(:,:,slice);
    imshow(rot90(tmp),[]); colormap('hot'); c = colorbar; c.Location = 'southoutside'; c.Color = [1 1 1];
end
h.Color = [0 0 0];
saveas(h,fullfile(outpath,'De_maps'),'fig');

close all

done = 'All the parametric maps have been saved into figures.';

disp(done);


end