function writeVideoToFile(Video, filename, frameRate)
%WRITEVIDEOTOFILE Writes video to file, using a specified frame rate
v = VideoWriter(strcat('Output/', filename));

if nargin == 3
    v.FrameRate = frameRate;
end

open(v);
writeVideo(v, Video);
close(v);

end

