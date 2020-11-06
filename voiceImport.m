
[voice1 , fs] = audioread('starter.wav');
[voice2 , fs] = audioread('embedadapt_1sample.wav');

lengthMin = min(length(voice1),length(voice2));
ds = 3;
frameNum = floor(lengthMin/frameLength / ds);
frameNumber = min(frameNum,frameNumber);
voice48K = [voice1(1:frameNumber*frameLength*ds) voice2(1:frameNumber*frameLength*ds)];
voice16K = voice48K(1:ds:end,:); % downsample from 48KS/s to 16 KS/s
voice16K = reshape(voice16K,frameLength,frameNumber,size(voice16K,2));
fs = fs/ds;

% sound(voice16K(:,1),fs)