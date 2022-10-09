













% Constants
cdir = "C:\Users\sfg99\Code\Summer Research\Matlab";
generated = "data_generated";


% Initial/Final Geometry (Generated)
g_files = cdir + "\" + generated + "\Xitwocompleted41";
final_iteration = 1;
g0_x = load(g_files + "X0.txt"); gF_x = load(g_files + "X" + final_iteration + ".txt");
g0_y = load(g_files + "Y0.txt"); gF_y = load(g_files + "Y" + final_iteration + ".txt");
g0_z = load(g_files + "Z0.txt"); gF_z = load(g_files + "Z" + final_iteration + ".txt");


% Plots
figure("Name", "Surface @ t = 0");
surf(g0_x,g0_y,g0_z, 'edgecolor', 'none');
% scatter3(g0_x,g0_y,g0_z);
axis equal
figure("Name", "Surface @ t = " + final_iteration);
surf(gF_x,gF_y,gF_z, 'edgecolor', 'none');
axis equal

figure("Name", "Contact Line @ t = " + final_iteration);  plot(gF_x(1,:),gF_y(1,:), 'k'); %.-k
hold on
plot(gF_x(end,:),gF_y(end,:), 'k');
plot(gF_x(:,1),gF_y(:,1), 'k');
plot(gF_x(:,end),gF_y(:,end), 'k');
hold off
axis equal

figure("Name", "Combine Graph");
surf(gF_x,gF_y,gF_z, 'edgecolor', 'none');
hold on
surf(g0_x,g0_y,g0_z, 'edgecolor', 'none', 'FaceColor','g', 'FaceAlpha',0.5);
hold off
axis equal