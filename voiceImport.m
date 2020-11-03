[voice1 , fs] = audioread('footsteps.wav');
[voice2 , fs] = audioread('heart.wav');
[voice3 , fs] = audioread('starter.wav');

time = 10; % second
nSample = time * fs;
voice48K = [voice1(1:nSample) voice2(1:nSample) voice3(1:nSample)];
ds = 3;
voice16K = voice48K(1:ds:end,:); % downsample from 48KS/s to 16 KS/s
fs = fs/ds;
voice16K(2,:) = circshift(voice16K(2,:),30);
% sound(voice16K(:,1),fs)