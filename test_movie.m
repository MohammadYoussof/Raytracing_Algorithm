x=1:256;
[x y] = meshgrid(x,x);
figure(1)
vidfile = VideoWriter('testmovie.mp4','MPEG-4');
open(vidfile);
for ind = 1:256
    z=sin(x*2*pi/ind)+cos(y*2*pi/ind);
   imagesc(z),colormap(hot) 
    drawnow
    F(ind) = getframe(gcf); 
    writeVideo(vidfile,F(ind));
end
close(vidfile)
