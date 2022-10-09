total_iterations = 24;

raw_dir = "C:\Users\sfg99\Code\Summer Research\Matlab\data_generated\";
ren_dir = "C:\Users\sfg99\Code\Summer Research\Matlab\render\";

for i = 2:total_iterations
   x = load(raw_dir + "Xitwocompleted41" + "X" + i + ".txt");
   y = load(raw_dir + "Xitwocompleted41" + "Y" + i + ".txt");
   z = load(raw_dir + "Xitwocompleted41" + "Z" + i + ".txt");
   fig = figure('visible', 'off');
   zlim([0 0.3]);
   axis equal;
   surf(x,y,z, 'edgecolor', 'none');
   saveas(fig, ren_dir + i, "png");
end