%--------------------------------------------------------------------------
%
% Example script for Signal_Generator
%
% Author        : dr.ir. Emanuel A.P. Habets
% Date          : 10-02-2008
%
% Related paper :
%
% Comment       :
%
%--------------------------------------------------------------------------

close all;
%clear;
set(0,'DefaultAxesFontSize',12)

%% Settings
c = 340;                         % Sound velocity in meters / second
LL = [4 5 6];                    % Room dimensions in meters ('golden' ratio)
beta = 0.2;                      % Reverberation time in seconds
nsample = 1024;                  % RIR length in samples
order = 2;                       % Reflection order
rp = [1.5 2.5-0.1 3; 1.5 2.5+0.1 3]; % Receiver positions in meter
cp = (rp(1,:) + rp(end,:))/2;
M = size(rp,1);
hop = 32;                        % Refresh rate of the AIR
sp = [2.5 0.5 3];                % Initial source position
type_of_movement = 'line';       % Source movement 'arc' or 'line'

%% Load anechoic sound source
[in, fs] = audioread('female_speech.wav');
if size(in,1) > size(in,2)
    in = in';
end
in = [in(1:4*fs) in(1:4*fs)];
len = length(in);

%% Generate source path
sp_path = zeros(len,3);
rp_path = zeros(len,3,M);

switch lower(type_of_movement)
    case 'arc'
        [theta phi r_d] = cart2sph(sp(1,1)-cp(1),sp(1,2)-cp(2),sp(1,3)-cp(3));
        theta = theta*180/pi;
        phi = phi*180/pi;
        
    case 'line'
        start_x = sp(1);
        start_y = sp(2);
        start_z = sp(3);
        stop_x = 2.5;
        stop_y = 4.5;
end

for ii = 1:hop:len
    switch lower(type_of_movement)
        % Calculate new source position (arc movement)
        case 'arc'
            [x_tmp, y_tmp, z_tmp] = sph2cart((theta+(ii*60)/len)*pi/180,phi*pi/180,r_d);
            sp = cp + [x_tmp y_tmp z_tmp];
            
            % Calculate new source position (line movement)
        case 'line'
            x_tmp = start_x + (ii*(stop_x-start_x)/len);
            y_tmp = start_y + (ii*(stop_y-start_y)/len);
            z_tmp = start_z;
            sp = [x_tmp y_tmp z_tmp];
    end
    
    % Store source path
    sp_path(ii:1:min(ii+hop-1,len),1) = sp(1);
    sp_path(ii:1:min(ii+hop-1,len),2) = sp(2);
    sp_path(ii:1:min(ii+hop-1,len),3) = sp(3);
    
    % Stationary receiver positions
    for mm=1:M
        rp_path(ii:1:min(ii+hop-1,len),1,mm) = rp(mm,1);
        rp_path(ii:1:min(ii+hop-1,len),2,mm) = rp(mm,2);
        rp_path(ii:1:min(ii+hop-1,len),3,mm) = rp(mm,3);    
    end
end

%% Generate microphone signal
[out,beta_hat] = signal_generator(in,c,fs,rp_path,sp_path,LL,beta,nsample,'o',order);

%% Plot source path and signals
figure(1);
plot3(rp(1,1),rp(1,2),rp(1,3),'x');
hold on;
for mm = 2:M
    plot3(rp(mm,1),rp(mm,2),rp(mm,3),'x');
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

%% Play generated microphone signal
% soundsc(out'./max(max(abs(out))),fs);