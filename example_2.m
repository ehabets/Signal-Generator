%--------------------------------------------------------------------------
%
% Example script for Signal_Generator
%
% Author        : dr.ir. Emanuel A.P. Habets
% Date          : 14-09-2011
%
% Related paper :
%
% Comment       :
%
%--------------------------------------------------------------------------

close all;
clear;
set(0,'DefaultAxesFontSize',12)

%% Settings
c = 340;           % Sound velocity in meters / second
LL = [4 5 6];      % Room dimensions in meters ('golden' ratio)
beta = 1;          % Reverberation time in seconds
order = 2;         % Reflection order
nsample = 256;     % RIR length in samples
M = 1;             % Number of microphones
sp = [2.5 4.5 3];  % Initial source position
hop = 32;          % Refresh rate of the AIR

% Source start and stop positions
start_x = sp(1);
start_y = sp(2);
stop_x = 2.5;
stop_y = 4.98;

%% Load anechoic sound source
[in, fs] = audioread('female_speech.wav');
%fs = 8000; in = randn(1,4*fs);
if size(in,1) > size(in,2)
    in = in';
end
in = [in(1:4*fs) in(1:4*fs)];
len = length(in);

%% Generate source and receiver paths
sp_path = zeros(len,3);
rp_path = zeros(len,3,M);

for ii = 1:hop:len
    % Calculate source position at sample ii
    x_tmp = start_x + (ii*(stop_x-start_x)/len);
    y_tmp = start_y + (ii*(stop_y-start_y)/len);
    
    % Store source path
    sp_path(ii:1:min(ii+hop-1,len),1) = x_tmp;
    sp_path(ii:1:min(ii+hop-1,len),2) = y_tmp;
    sp_path(ii:1:min(ii+hop-1,len),3) = sp(3);
    
    % Store receiver path (fixed offset)
    for mm = 1:M
        rp_path(ii:1:min(ii+hop-1,len),1,mm) = x_tmp-mm*0.5;
        rp_path(ii:1:min(ii+hop-1,len),2,mm) = y_tmp;
        rp_path(ii:1:min(ii+hop-1,len),3,mm) = sp(3);
    end    
end

%% Generate microphone signal
[out,beta_hat] = signal_generator(in,c,fs,rp_path,sp_path,LL,beta,nsample,'o',order);

%% Generate start and end RIR
% [h1,beta_hat] = rir_generator(c,fs,rp_path(1,:,1),sp_path(1,:,1),LL,beta,nsample,'o',order);
% [h2,beta_hat] = rir_generator(c,fs,rp_path(end,:,1),sp_path(end,:,1),LL,beta,nsample,'o',order);
% figure(10);plot([h1 ; h2].');

%% Plot source path and signals
figure(1);
plot3(rp_path(:,1,1),rp_path(:,2,1),rp_path(:,3,1),'x');
hold on;
for mm = 2:M
    plot3(rp_path(:,1,mm),rp_path(:,2,mm),rp_path(:,3,mm),'x');
end
plot3(sp_path(:,1),sp_path(:,2),sp_path(:,3),'r.');
axis([0 LL(1) 0 LL(2) 0 LL(3)]);
grid on;
box on;
axis square;
hold off;

figure(2)
t = 0:1/fs:(length(in)-1)/fs;
subplot(211); plot(t,in); title('in(n)');xlabel('Time [Seconds]');ylabel('Amplitude');
subplot(212); plot(t,out'); title('out(n)');xlabel('Time [Seconds]');ylabel('Amplitude');
grid on;

%% Play generated microphone signal
% soundsc(out'./max(max(abs(out))),fs);