function writeVideoToFile(Video, filename, frameRate)
%WRITEVIDEOTOFILE Summary of this function goes here
%   Detailed explanation goes here

v = VideoWriter(strcat('Output/', filename));

if nargin == 3
    v.FrameRate = frameRate;
end

open(v);
writeVideo(v, Video);
close(v);

end

