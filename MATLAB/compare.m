% Constants
cdir = "C:\Users\sfg99\Code\Summer Research\Matlab";
source = "data_previous"; % data_previous, data_original
generated = "data_generated";


% Initial/Final Geometry (Source)
e_files = cdir + "\" + source + "\Xitwocompleted41";
e0_x = load(e_files + "X0.txt"); eF_x = load(e_files + "X1.txt");
e0_y = load(e_files + "Y0.txt"); eF_y = load(e_files + "Y1.txt");
e0_z = load(e_files + "Z0.txt"); eF_z = load(e_files + "Z1.txt");

% Initial/Final Geometry (Generated)
g_files = cdir + "\" + generated + "\Xitwocompleted41";
g0_x = load(g_files + "X0.txt"); gF_x = load(g_files + "X1.txt");
g0_y = load(g_files + "Y0.txt"); gF_y = load(g_files + "Y1.txt");
g0_z = load(g_files + "Z0.txt"); gF_z = load(g_files + "Z1.txt");

% Initial/Final Geometry Difference
d0_x = g0_x - e0_x; dF_x = gF_x - eF_x;
d0_y = g0_y - e0_y; dF_y = gF_y - eF_y;
d0_z = g0_z - e0_z; dF_z = gF_z - eF_z;


% Statistics Calculations
N = numel(dF_x);
mean_diff_x = sum(dF_x, "all") / N;
mean_diff_y = sum(dF_y, "all") / N;
mean_diff_z = sum(dF_z, "all") / N;
mean_absd_x = sum(abs(dF_x), "all") / N;
mean_absd_y = sum(abs(dF_y), "all") / N;
mean_absd_z = sum(abs(dF_z), "all") / N;
% strd_dev_x = sqrt(sum( (dF_x - mean_diff_x * ones(21,21)) .* (dF_x - mean_diff_x * ones(21,21)), "all")) / sqrt(N + 1);
% strd_dev_y = sqrt(sum( (dF_y - mean_diff_y * ones(21,21)) .* (dF_y - mean_diff_y * ones(21,21)), "all")) / sqrt(N + 1);
% strd_dev_z = sqrt(sum( (dF_z - mean_diff_z * ones(21,21)) .* (dF_z - mean_diff_z * ones(21,21)), "all")) / sqrt(N + 1);
disp("Mean difference     (" + mean_diff_x + ", " + mean_diff_y + ", " + mean_diff_z + ")")
disp("Abs Mean difference (" + mean_absd_x + ", " + mean_absd_y + ", " + mean_absd_z + ")")
% disp("Standard Deviation  (" + strd_dev_x  + ", " + strd_dev_y  + ", " + strd_dev_z  + ")")


% Plots

% figure("Name", "Expected @ t = 0");      surf(e0_x,e0_y,e0_z, 'edgecolor', 'none'); axis equal;
% figure("Name", "Expected @ t = 1200");   surf(eF_x,eF_y,eF_z, 'edgecolor', 'none'); axis equal;
figure("Name", "Generated @ t = 0");     surf(g0_x,g0_y,g0_z, 'edgecolor', 'none'); axis equal;
figure("Name", "Generated @ t = 1200");  surf(gF_x,gF_y,gF_z, 'edgecolor', 'none'); axis equal;
figure("Name", "Difference @ t = 0");    surf(d0_x,d0_y,d0_z, 'edgecolor', 'none'); axis equal;
figure("Name", "Difference @ t = 1200"); surf(dF_x,dF_y,dF_z, 'edgecolor', 'none'); axis equal;

% Misc comments
% scatter3 for 3D scatter plot

% figure("Name", "Double Plot @ t = 1200")
% surf(gF_x,gF_y,gF_z, 'edgecolor', 'none', 'FaceColor','g', 'FaceAlpha',0.25)
% hold on
% surf(eF_x,eF_y,eF_z, 'edgecolor', 'none', 'FaceColor','r', 'FaceAlpha',0.75)
% hold off
% axis equal

% quiver3;
% xlabel('x')