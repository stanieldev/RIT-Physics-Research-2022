total_iterations = 24;


ren_dir = "C:\Users\sfg99\Code\Summer Research\Matlab\render\";

writerObj = VideoWriter('animation.avi');
writerObj.FrameRate = 2;

open(writerObj);
for i=2:total_iterations
    I = imread(ren_dir + i + ".png") ; 
    imshow(I) ; 
    F = getframe(gcf);
    writeVideo(writerObj, F);
end
close(writerObj);