function playSound(type)
%PLAYSOUND Summary of this function goes here
%   Detailed explanation goes here

if type == "Complete" || type == "complete"
    WarnWave = [sin(1:.6:400), sin(1:.7:400), sin(1:.4:400)];
    Audio = audioplayer(WarnWave, 22050);
    play(Audio);
end

end

