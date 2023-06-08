clc
clear
close all
%% Importing data

% path to .mat file with images
path_im_mac = 'sim_data_init_plus_final.mat';
path_im_windows = 'sim_data_init_plus_final.mat';

% Path to .mat file with I vs q graphs
path_Ivsq_mac = 'I_vs_q_lab_sims.mat';
path_Ivsq_windows = 'I_vs_q_lab_sims.mat';

% Path to .mat file for boundry data
path_bound_mac = 'sim_results_M_s_boundry_lab.mat';
path_bound_windows = 'sim_results_M_s_boundry_lab.mat';

% Path to xlsx file with shift values from 2D fits
path_shift_data_mac = 'velocity_calc.xlsx';
path_shift_data_windows = 'velocity_calc.xlsx';

% Path to export xlsx file
path_export_data_mac = '/Exported_results/';
path_export_data_windows = '\Exported_results\';

% Setting directory depending on the OS
if ismac
    path_im = path_im_mac;
    path_Ivsq = path_Ivsq_mac;
    path_bound = path_bound_mac;
    path_shift_data = path_shift_data_mac;
    path_export_data = path_export_data_mac;
elseif isunix
    path_im = path_im_linux;
    path_Ivsq = path_Ivsq_linux;
    path_bound = path_bound_linux;
    path_shift_data = path_shift_data_linux;
    path_export_data = path_export_data_linux;
elseif ispc
    path_im = path_im_windows;
    path_Ivsq = path_Ivsq_windows;
    path_bound = path_bound_windows;
    path_shift_data = path_shift_data_windows;
    path_export_data = path_export_data_windows;
else
    disp('Platform not supported')
end

% Loading image data from micromag simulations
all_im_data = load(path_im, 'save_mat');
all_im_data = cell2mat(struct2cell(all_im_data));

% Loading Ivsq data from micromag simulations
all_Ivsq_data = load(path_Ivsq, 'q_save_data');
all_Ivsq_data = cell2mat(struct2cell(all_Ivsq_data));

% Loading boundry data from micromag simulations
load(path_bound, 'B_Lab');
all_bound_data = B_Lab;

%% Setting parameters
pix_size = 2; % pix size in nm
image_dim = size(all_im_data(:, :, 3)); % image size in matrix [n x m]
area_threshold_num = 10; % Threshold for min area of a feature in diff image
feret_minor_threshold_num_high = 35; % Upper threshold for feret minor dimension
feret_minor_threshold_num_low = 2; % Lower threshold for feret minor dimension
q_threshold = 0.05; % Setting upper limit for Ivsq data x (q) axis
%% Seperating data
% Creating seperate x, y, z variables (for image) for easy handling
X = all_im_data(:, :, 1); % X axis of images
Y = all_im_data(:, :, 2); % Y axis of images
im_inital = all_im_data(:, :, 3); % Initail image ()
im_final = all_im_data(:, :, 4); % Final image after 

% Creating seperate x, y, z variables (for image) for easy handling
q = all_Ivsq_data(all_Ivsq_data(:, 1)<=q_threshold, 1); % q for Ivsq
I_initial = all_Ivsq_data(all_Ivsq_data(:, 1)<=q_threshold, 2); % I initial
I_final = all_Ivsq_data(all_Ivsq_data(:, 1)<=q_threshold, 3); % I final post sim evolution

%% Taking an FFT
fft_resolution_X = 2*pi/(size(X, 1)*2)*2;
fft_resolution_Y = 2*pi/(size(Y, 1)*2)*2;

X_fft = linspace(-fft_resolution_X*size(X, 1)/2, fft_resolution_X*size(X, 1)/2, size(X, 1));
Y_fft = linspace(-fft_resolution_Y*size(Y, 1)/2, fft_resolution_Y*size(Y, 1)/2, size(Y, 1));

[X_fft, Y_fft] = meshgrid(X_fft, Y_fft);
X_fft = X_fft';
Y_fft = Y_fft';

figure
pcolor(X, Y, im_inital);
shading flat;
colormap gray;
axis square;
ax=gca;ax.FontSize=16;
xlabel('X (nm)');
ylabel('Y (nm)');
clim([-1 1]);
colorbar;
title('Initial image')



fft_one_img = (fft2(im_inital));
fft_2D_img = abs(fftshift(fft_one_img));
% figure
% imagesc(test(424:600, 424:600))
% axis equal
% axis square

figure
pcolor(X_fft, Y_fft, fft_2D_img);
shading flat;
colormap jet;
axis square;
ax=gca;ax.FontSize=16;
xlabel('X (1/nm)');
ylabel('Y (1/nm)');
xlim([-0.15, 0.15]);
ylim([-0.15, 0.15]);
% clim([-1 1]);
colorbar;
title('FFT image')

%%
fix_axis_value = 512;

figure;
% imagesc(fft_2D_img(480:545, 480:545));
imagesc(fft_2D_img);
axis equal
axis square

x_line = [512 512+25];
y_line = [fix_axis_value fix_axis_value];

figure
c = improfile(fft_2D_img, x_line, y_line);
x_axis = linspace(0, size(c, 1)-1, size(c, 1))*fft_resolution_X;
plot(x_axis, c);
xlabel('X (1/nm)');
ylabel('Intensity (arb)');
title('X+')

x_line = [512-25 512];
y_line = [fix_axis_value fix_axis_value];

figure
c = improfile(fft_2D_img, x_line, y_line);
x_axis = linspace(0, size(c, 1)-1, size(c, 1))*fft_resolution_X;
plot(x_axis, c);
xlabel('X (1/nm)');
ylabel('Intensity (arb)');
title('X-')

x_line = [fix_axis_value fix_axis_value];
y_line = [512 512+25];

figure
c = improfile(fft_2D_img, x_line, y_line);
x_axis = linspace(0, size(c, 1)-1, size(c, 1))*fft_resolution_X;
plot(x_axis, c);
xlabel('X (1/nm)');
ylabel('Intensity (arb)');
title('Y+')

x_line = [fix_axis_value fix_axis_value];
y_line = [512-25 512];

figure
c = improfile(fft_2D_img, x_line, y_line);
x_axis = linspace(0, size(c, 1)-1, size(c, 1))*fft_resolution_X;
plot(x_axis, c);
xlabel('X (1/nm)');
ylabel('Intensity (arb)');
title('Y-')
%% Looking at boundry data and raw image

figure;
pcolor(X,Y,im_final);shading flat;colormap gray;
hold on;
for k = 1:length(B_Lab)
    boundary = B_Lab{k};
    curve_for_line = LineCurvature2D(boundary*2);
    plot(boundary(:,1)*2,boundary(:,2)*2,'Color','red','LineWidth',1.5);
%     scatter(boundary(:,1)*2,boundary(:,2)*2, 1, curve_for_line);
end
colorbar
hold off;
axis square;
ax=gca;ax.FontSize=16;
xlabel('X (nm)');ylabel('Y (nm)');
xlim([100 1000])
ylim([1000 1900])
zlabel('M_Z (arb)')
clim([-1 1]);%colorbar;

% figure;
% pcolor(X,Y,im_inital);shading flat;colormap gray;
% hold on;
% for k = 1:length(B_Lab)
%     boundary = B_Lab{k};
%     curve_for_line = LineCurvature2D(boundary*2);
%     plot(boundary(:,1)*2,boundary(:,2)*2,'Color','red','LineWidth',1.5);
% %     scatter(boundary(:,1)*2,boundary(:,2)*2, 1, curve_for_line);
% end
% colorbar
% hold off;
% axis square;
% ax=gca;ax.FontSize=16;
% xlabel('X (nm)');ylabel('Y (nm)');
% caxis([-1 1]);%colorbar;

% % Trying 2D FFT
% Y = fft2(im_inital);
% figure
% imagesc(abs(fftshift(Y)))
%% Calculating FFT q change in simulations using gauss fits
% fitting initial and final I vs q data with Gaussians
fit_initial = fit(q, I_initial, 'gauss1');
fit_final = fit(q, I_final, 'gauss1');

% fit_initial.b1
% fit_final.b1

fprintf('q initial from FFT: %s 1/nm \n', num2str(round(fit_initial.b1, 5)));
fprintf('q final from FFT: %s 1/nm \n', num2str(round(fit_final.b1, 5)));

% Calculating mean shift
mean_FFT_shift = pi*(fit_initial.b1-fit_final.b1)/(fit_initial.b1*fit_final.b1);
% disp('Change in domain size:', sim_FFT_delta_domain)
fprintf('Change in domain size from FFT: %s nm \n', num2str(round(mean_FFT_shift, 3)));

fract_FFT_shift = 2*(fit_initial.b1-fit_final.b1)/(fit_initial.b1+fit_final.b1);
% disp('Change in domain size:', sim_FFT_delta_domain)
fprintf('Frac shift in q from FFT: %s 1/nm \n', num2str(round(fract_FFT_shift, 5)));

% Plotting gauss fits on FFT for sim
figure
hold on
box on
plot(fit_initial,'b-', q, I_initial, 'b*');
xline(fit_initial.b1, 'b', 'HandleVisibility','off')
plot(fit_final, q, I_final, 'rx');
xline(fit_final.b1, 'r', 'HandleVisibility','off')
xlabel('q (1/nm)');
ylabel('Int Intensity (arb)');
legend('Data\_equlibrium', 'Fit\_equlibrium', 'Data\_200\_ps', 'Fit\_200\_ps')
hold off

% P0 = [1000 0.03 0.01 1000];
% % BOUNDS = [LB1 LB2 LB3 LB4; UB1 UB2 UB3 UB4];
% 
% [yprime2 params2 resnorm2 residual2] = lorentzfit(q,I_final,[],[],'3c');
% figure; plot(q,I_final,'b.','LineWidth',2)
% hold on; plot(q,yprime2,'r-','LineWidth',2)
%% Plotting initial image
% figure
% pcolor(X,Y,im_inital);shading flat;colormap gray;
% axis square;
% ax=gca;ax.FontSize=16;
% xlabel('X (nm)');
% ylabel('Y (nm)');
% clim([-1 1]);
% colorbar;
% 
% %% Plotting final image
% figure
% pcolor(X, Y, im_final);shading flat;colormap gray;
% axis square;
% ax=gca;ax.FontSize=16;
% xlabel('X (nm)');
% ylabel('Y (nm)');
% clim([-1 1]);
% colorbar;

%% Binerize data and plotting
% close all

% Bineriseing initial and final image
im_inital_bin = binerise_data(im_inital);
im_final_bin = binerise_data(im_final);
% Calculating the difference
im_diff = im_final_bin - im_inital_bin;

% Plotting initial
figure
pcolor(X, Y, im_inital_bin);
shading flat;
colormap gray;
axis square;
ax=gca;ax.FontSize=16;
xlabel('X (nm)');
ylabel('Y (nm)');
clim([0 1]);
colorbar;
title('Initial image')

% Plotting final
figure
pcolor(X, Y, im_final_bin);
shading flat;
colormap gray;
axis square;
ax=gca;ax.FontSize=16;
xlabel('X (nm)');
ylabel('Y (nm)');
clim([0 1]);
colorbar;
title('Final image')

% Plotting difference
figure
pcolor(X, Y, im_diff);
shading flat;
colormap gray;
axis square;
ax=gca;ax.FontSize=16;
xlabel('X (nm)');
ylabel('Y (nm)');
clim([-1 1]);
colorbar;
title('Diff image')

%% Using regionprops to analyse the data greater than 0
close all

% Only looking at the + regions in diff image
im_diff_p = im_diff>0;

% figure
% pcolor(X, Y, im_diff_p.*1);
% shading flat;
% colormap gray;
% axis square;
% ax=gca;ax.FontSize=16;
% xlabel('X (nm)');
% ylabel('Y (nm)');
% clim([0 1]);
% colorbar;
% title('Final image')

% Using regionprops to get the properties of each structure
stats_p = regionprops(im_diff_p, 'all');
% Getting stats from reginprops
area_p = cat(1, stats_p.Area);
centroids = cat(1,stats_p.Centroid);
minor_axis = cat(1, stats_p.MinorAxisLength);
major_axis_p = cat(1, stats_p.MajorAxisLength);

% Uing feret minor axis for threshold. Not using anymore for new algo
feret_minor = cat(1, stats_p.MinFeretDiameter);
feret_minor_threshold_high = (feret_minor<feret_minor_threshold_num_high);
feret_minor_threshold_low = (feret_minor>feret_minor_threshold_num_low);

% Mask areas below area threhold
area_threshold = (area_p>area_threshold_num);
final_mask_p = area_threshold;% & feret_minor_threshold_high & feret_minor_threshold_low;

% Showing centroids on the image
% figure
% imshow(im_diff_p)
% hold on
% plot(centroids(final_mask_p,1),centroids(final_mask_p,2),'b*')
% hold off

data_width_p = minor_axis(final_mask_p);
stats_for_p_shift = datastats(data_width_p);

%% Plotting the ellipes on image for positive shift values in diff image
close all
figure
imshow(im_diff_p)
t = linspace(0,2*pi,50);
hold on
for k = 1:length(stats_p)
    if final_mask_p(k) == true
        a = stats_p(k).MajorAxisLength/2;
        b = stats_p(k).MinorAxisLength/2;
        Xc = stats_p(k).Centroid(1);
        Yc = stats_p(k).Centroid(2);
        phi = deg2rad(-stats_p(k).Orientation);
        x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
        y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
        plot(x,y,'r','Linewidth',2)
    end
end
hold off

%% Using regionprops to analyse the data less than 0
% close all
im_diff_n = im_diff<0;

% figure
% pcolor(X, Y, im_diff_p.*1);
% shading flat;
% colormap gray;
% axis square;
% ax=gca;ax.FontSize=16;
% xlabel('X (nm)');
% ylabel('Y (nm)');
% clim([0 1]);
% colorbar;
% title('Final image')

% Using regionprops to get the properties of each structure
stats_n = regionprops(im_diff_n, 'all');
area_n = cat(1, stats_n.Area);
centroids = cat(1,stats_n.Centroid);
minor_axis = cat(1, stats_n.MinorAxisLength);
major_axis_n = cat(1, stats_n.MajorAxisLength);

feret_minor = cat(1, stats_n.MinFeretDiameter);
feret_minor_threshold_high = (feret_minor<feret_minor_threshold_num_high);
feret_minor_threshold_low = (feret_minor>feret_minor_threshold_num_low);

area_threshold = (area_n>area_threshold_num);
final_mask_n = area_threshold;% & feret_minor_threshold_high & feret_minor_threshold_low;
figure
imshow(im_diff_n)
hold on
plot(centroids(final_mask_n,1),centroids(final_mask_n,2),'b*')
hold off

data_width_n = minor_axis(final_mask_n);
stats_for_n_shift = datastats(data_width_n);

%% Plotting the ellipes on image for n
% close all
% figure
% imshow(im_diff_n)
% t = linspace(0,2*pi,50);
% hold on
% for k = 1:length(stats_n)
%     if final_mask_n(k) == true
%         a = stats_n(k).MajorAxisLength/2;
%         b = stats_n(k).MinorAxisLength/2;
%         Xc = stats_n(k).Centroid(1);
%         Yc = stats_n(k).Centroid(2);
%         phi = deg2rad(-stats_n(k).Orientation);
%         x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
%         y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
%         plot(x,y,'r','Linewidth',2)
%     end
% end
% hold off

%% Plotting the ellipes on image for whole image
close all

figure
pcolor(im_diff);
shading flat;
colormap gray;
axis square;
ax=gca;ax.FontSize=16;
xlabel('X (nm)');
ylabel('Y (nm)');
clim([-1 1]);
colorbar;
title('Diff image')
% imshow(im_diff)

t = linspace(0,2*pi,100);
extra_wiggle = 4;
hold on
for k = 1:length(stats_p)
    if final_mask_p(k) == true
        a = stats_p(k).MajorAxisLength/2;
        b = stats_p(k).MinorAxisLength/2;
        Xc = stats_p(k).Centroid(1);
        Yc = stats_p(k).Centroid(2);
        phi = deg2rad(-stats_p(k).Orientation);
        x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
        y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
        plot(x,y,'r','Linewidth',2)
        minor_angle = deg2rad(-stats_p(k).Orientation+90);
        m = tan(minor_angle);
        c = Yc-m*Xc;
        x_range = [(Xc-(stats_p(k).MinFeretDiameter/2+extra_wiggle)*cos(minor_angle)), Xc+(stats_p(k).MinFeretDiameter/2+extra_wiggle)*cos(minor_angle)];
        x = linspace(x_range(1), x_range(2), 100);
        y = m*x+c;
%         plot(Xc, Yc, 'g*')
        plot(x, y, 'r', 'LineWidth',2)
        listcut_p{k, 1} = improfile((im_diff_p), x, y);
        line_len = sqrt(((max(x)-min(x))*pix_size)^2 + ((max(y)-min(y))*pix_size)^2);
        norm_width = sum(listcut_p{k, 1},'all', 'omitnan')/size(listcut_p{k, 1}, 1);
        % This is the final width
        final_width_p(k, 1) = line_len*norm_width;
        % List of areas
        final_area_p(k, 1) = area_p(k, 1);
    else
        listcut_p{k, 1} = nan;
        final_width_p(k, 1) = nan;
        final_area_p(k, 1) = nan;
    end
end

for k = 1:length(stats_n)
    if final_mask_n(k) == true
        a = stats_n(k).MajorAxisLength/2;
        b = stats_n(k).MinorAxisLength/2;
        Xc = stats_n(k).Centroid(1);
        Yc = stats_n(k).Centroid(2);
        phi = deg2rad(-stats_n(k).Orientation);
        x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
        y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
        plot(x,y,'b','Linewidth',2)
        minor_angle = deg2rad(-stats_n(k).Orientation+90);
        m = tan(minor_angle);
        c = Yc-m*Xc;
        x_range = [(Xc-(stats_n(k).MinFeretDiameter/2+extra_wiggle)*cos(minor_angle)), Xc+(stats_n(k).MinFeretDiameter/2+extra_wiggle)*cos(minor_angle)];
        x = linspace(x_range(1), x_range(2), 100);
        y = m*x+c;
%         plot(Xc, Yc, 'g*')
        plot(x, y, 'b', 'LineWidth',2)
        listcut_n{k, 1} = improfile((im_diff_n), x, y);
        line_len = sqrt(((max(x)-min(x))*pix_size)^2 + ((max(y)-min(y))*pix_size)^2);
        norm_width = sum(listcut_n{k, 1},'all', 'omitnan')/size(listcut_n{k, 1}, 1);
        % This is the final width
        final_width_n(k, 1) = line_len*norm_width;
        % List of areas
        final_area_n(k, 1) = area_n(k, 1);
    else
        listcut_n{k, 1} = nan;
        final_width_n(k, 1) = nan;
        final_area_n(k, 1) = nan;
    end
end

% for k = 1:length(stats_n)
%     if final_mask_n(k) == true
%         a = stats_n(k).MajorAxisLength/2;
%         b = stats_n(k).MinorAxisLength/2;
%         Xc = stats_n(k).Centroid(1);
%         Yc = stats_n(k).Centroid(2);
%         phi = deg2rad(-stats_n(k).Orientation);
%         x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
%         y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
%         plot(x,y,'r','Linewidth',2)
%     end
% end
hold off


%% Final calculations
mean_area_shift_p = rms(final_area_p*pix_size^2, 'all', 'omitnan');
mean_area_shift_n = rms(final_area_n*pix_size^2, 'all', 'omitnan');

test = [final_area_n*pix_size^2;final_area_p*pix_size^2];

mean_area_shift = rms(test, 'all', 'omitnan');
fprintf('Change in area from Image: %s nm^2 \n', num2str(round(mean_area_shift, 3)));

test = sum(final_area_p, 'all', 'omitnan')+sum(final_area_n, 'all', 'omitnan')

frac_area_shift = (sum(final_area_p, 'all', 'omitnan')+sum(final_area_n, 'all', 'omitnan'))/(prod(image_dim));
fprintf('Frac change in area from Image: %s arb \n', num2str(round(frac_area_shift, 5)));

mean_im_shift_p = rms(final_width_p, 'all', 'omitnan');
mean_im_shift_n = rms(final_width_n, 'all', 'omitnan');

mean_im_shift = rms([final_width_p;final_width_n], 'all', 'omitnan');
fprintf('Change in domain size from Image: %s nm\n', num2str(round(mean_im_shift, 3)));

curvature_density = (sum(final_mask_p, 'all') + sum(final_mask_n, 'all'))/(prod(image_dim)*pix_size^2/(1000^2));
fprintf('Curvature density for sim is: %s um^-2\n', num2str(round(curvature_density, 2)));

total_wall_len = (sum(major_axis_p*pix_size, 'all') + sum(major_axis_n*pix_size, 'all'));
fprintf('Wall len: %s um\n', num2str(round(total_wall_len/1000, 2)));

frac_shift_from_mean = mean_area_shift*curvature_density/1000000

multiplier = frac_shift_from_mean/fract_FFT_shift

multiplier_2 = frac_area_shift/fract_FFT_shift*curvature_density/1000000

multiplier_3 = mean_area_shift/(prod(image_dim*pix_size))/fract_FFT_shift
%% Final calc
FERMI_experiment_data = readtable(path_shift_data);
Fluence_list = FERMI_experiment_data.Fluence_mJ_cm2_;

domain_size_change_scattering_experiment_list = FERMI_experiment_data.dSize_nm_;

domain_size_change_scattering = domain_size_change_scattering_experiment_list;%5.4;
curvature_density_MFM = 60.76; % Curvature density in #/um^2

final_velocity = (domain_size_change_scattering)*(mean_im_shift/mean_FFT_shift)...
    *(curvature_density/curvature_density_MFM)/0.25;

fprintf('Domain wall velocity is: %s km/s\n', num2str(round(final_velocity(end), 2)));

%% Plotting fluence vs velocity
% close all
figure
plot(Fluence_list, final_velocity,'-*', LineWidth=2)
box on
ax = gca;
ax.FontSize = 14;
xlabel('Fluence (mJ/cm^2)', fontsize = 16)
ylabel('Velocity (km/s)', fontsize = 16)
%% Finding curvature and wall length
% close all
% im_initial_bin_upscale = imresize(im_inital_bin, 2);
% im_edge_initial_sig_low = edge(im_initial_bin_upscale,'Canny',[], 0.5);
% im_edge_initial_sig_high = edge(im_initial_bin_upscale,'Canny',[], 3);
% figure
% imshowpair(im_initial_bin_upscale, im_edge_initial_sig_low,'montage')
% figure
% imshowpair(im_edge_initial_sig_low, im_edge_initial_sig_high,'montage')
table_var_names = {'Ms_change (%)', 'q_i_fft (1/nm)', 'q_f_fft (1/nm)', 'dA_area_RMS_im (nm^2)', 'frac_area_change (arb)','d_shift_RMS_im (nm)', 'curvature_density (um^-2)', 'Avg_frac_shift_FFT', 'Multip'};

data_export_cell = {40, fit_initial.b1, fit_final.b1, mean_area_shift, frac_shift_from_mean, mean_im_shift, curvature_density, fract_FFT_shift, multiplier};

export_table = cell2table(data_export_cell, 'VariableNames', table_var_names);

filename = strcat(path_export_data, 'velocity_calc_data_RJ_final.xlsx');
% writetable(export_table,filename,'Sheet',1)
%% Binerizing data function

function [binerized_data] = binerise_data(mat)

mask_p = (mat>=0).*1;
mask_n = (mat<0).*0;
binerized_data = mask_p + mask_n;
end