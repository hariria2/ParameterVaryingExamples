function SaveVid(frames, filename)

writerObj = VideoWriter(filename);
open(writerObj);

for k = 1:length(frames) 
   writeVideo(writerObj,frames(k));
end

close(writerObj);