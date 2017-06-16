


outputVideo = VideoWriter(fullfile('video','ais.avi'));
open(outputVideo)

for ii = 1:length(imageNames)
   img = imread(sprintf('video/image%d.jpg', ii));
   writeVideo(outputVideo,img)
end

close(outputVideo)




