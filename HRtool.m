%% HR tool for ECG and POX
% by Javier Jaimovich (2012)
% version 3.4 (updated 21/06/2012)
%
% [HR, meanHR, Q] = HRtool(ECG, SR, debug, par, input)
%
% Returns a matrix with HR (time, bpm) and average HR
% Q is the level of confidence of the extraction (%)
%
% ECG: Input signal (ECG or POX)
% SR: Sampling Rate
% debug: 0 for no debugging, 1 to debug (includes plot)
% par is a vector with [minHR maxHR change_ratio R]
% being:
%   minHR: minimum HR limit
%   maxHR: maximum HR limit
%   R: Ratio between beat and std of ECG for beat threshold (only for ECG)
%   change_ratio: Maximum change ratio between beats (e.g. 40% -> 0.4)
%   
%   Default values:
%       debug = 0
%       minHR = 50
%       maxHR = 140
%       change_ratio = 0.3 (30%)
%       R = 3
%       input = 'ECG'
%
%   note: For POX signals HRtool will use 1st sample above threshold for
%   calculations. With ECG signals, HRtool will use max absolute value of
%   each waveform beat.
%   
%   ***ISSUES***
%   - Q is not considering signals that stop before the end of the file
%   (e.g. healthy pox signal for the first 30 sec of a 90 sec file). The
%   calculation should be: find de time of the last pulse and if it's more
%   than minHR distance from the end of the file, deduct points to Q
%   proportionally to that distance (e.g. minHR = 60, last heart beat at 5
%   seconds from the end of a 90 sec file, then Q = Q - 100*(5-1)/30; Look
%   at eimCorrectHR_Q.m file


function [HR, meanHR, Q] = HRtool(varargin)
%% Set-Up
version = 'v3.3';

%Definitions & Default values
startBPM = 70; %First beat in signal
T_max = 5; %Maximum number of replaced beats before resetting previous HR
n_debug = 0; %default debug
min_HR = 50;
max_HR = 140;
change_ratio = 0.3;
R = 3;
signal = 'ECG';
CAL_LENGTH = 10; %Calibration length in seconds

%% Check input parameters

switch nargin
    case 2
        ECG = varargin{1};
        SR = varargin{2};
    case 3
        ECG = varargin{1};
        SR = varargin{2};
        n_debug = varargin{3};
    case 4
        ECG = varargin{1};
        SR = varargin{2};
        n_debug = varargin{3};
        par = varargin{4};
        if length(par) ~= 4;
            error('There should be 4 parameters in vector'); end
        min_HR = par(1);
        max_HR = par(2);
        change_ratio = par(3);
        R = par(4);
    case 5
        ECG = varargin{1};
        SR = varargin{2};
        n_debug = varargin{3};
        par = varargin{4};
        if length(par) ~= 4;
            error('There should be 4 parameters in vector'); end
        min_HR = par(1);
        max_HR = par(2);
        change_ratio = par(3);
        R = par(4);
        signal = varargin{5};
        if (strcmpi(signal,'ECG')|strcmpi(signal,'POX'))~=1;
            error('input option should be "ECG" or "POX"'); end
    otherwise
        error('Check input parameters');
end

if SR <= 10; error('SR must be greater than 10'); end
if n_debug ~= 0 && n_debug ~= 1
    error('debug must be 0 or 1')
end

N = length(ECG);
if (n_debug==1);disp(['HRcalc ' version]);end

%Check that signal is longer than calibration length
if (N/SR) < CAL_LENGTH, error('Signal must be longer than %d (s)',CAL_LENGTH); end

%% Pre-Processing
T = 1/SR; %Interval

N = length(ECG); %Lenght of input signal
t = (0:N-1).*T; %create timeline

ECG_D = detrend(ECG); %Detrend ECG signal

switch strcmpi(signal,'ECG')
    case 1 %for ECG signal
        %High Pass filter (Kaiser Window)
        Fc = [2 3];      %Cutoff Frequencies (Hz)
        
        %Filter Coefficients
        [n,Wn,beta,ftype] = kaiserord(Fc,[0 1],[0.05 0.01],SR); %Kaiser
        Num_HR = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale'); %Create coefficients
        
        filt_ECG = filter(Num_HR,1,ECG_D); %Filter
        %Compensate filter latency
        left = floor((length(Num_HR)-1)/-2);
        ECG_A = abs(circshift(filt_ECG, left)); %Get absolute value
        
    case 0 %for POX signal
        thr = std(ECG_D)*1.2; %STD + 20%
        ECG_A = max(thr,ECG_D)-thr; %only values over threshold
        ECG_A(find(ECG_A)) = 1; %force all values ~= 0 to be ones
end
        
ECG_S = zeros(1,N); %preallocate memory for ECG_S: Values over threshold

%% Find samples above threshold level

b = floor(2.*SR); %Buffer Size of 2 seconds

%Filter values that are over threshold
if strcmpi(signal,'ECG')==1
    for i = 1:b:N
        %Check if buffer size is greater than signal at the end of analysis
        if (N - i) < b
            L = N-i;
        else
            L = b; %L is the actual buffer size used
        end
        %Z-normalization of detrended buffer
        ECG_Z(i:i+L) = zscore(ECG_A(i:i+L));
        %Set threshold for buffer (R standard deviations over mean)
        S(i:i+L) = mean(ECG_A(i:i+L)) + R.*std(ECG_A(i:i+L));
        for m = i:i+L %Extract values over threshold
            if ECG_A(m) > S(i);
                ECG_S(m) = ECG_Z(m);
            else
                ECG_S(m) = 0;
            end
        end
    end
else
    ECG_S = ECG_A;
    S(1:length(ECG_S)) = 0; %POX threshold
end

%% 1st Beat Calibration

%Get estimate of mean HR during 1st beats for setting initial BPM
CAL_N = floor(CAL_LENGTH*SR); %calibration length in samples
cal_buff = ECG_A(1:CAL_N); % calibration buffer
cal_thr = max(cal_buff)*0.6; % 60% of max value
cal_buff = max(cal_buff,cal_thr)-cal_thr; %filter values above thr
i=1; j=1; cal_b=zeros(1,1); %loop variables

while i<CAL_N %search for values above threshold
    if cal_buff(i)>0
        cal_b(j)=i; j=j+1; i=i+floor(SR*60/(max_HR));
    else i=i+1;
    end
end

% Calculate initial BPM and check if it is within ranges
if length(cal_b)>1;
    tempBPM = SR*60./diff(cal_b); %vector of all HR found in buffer
    tempBPM = sort(tempBPM); %sorted in ascending order 
    tempBPM = mean(tempBPM(max(floor(length(tempBPM)*.25),1):...
        ceil(length(tempBPM)*.75))); %extract mean from middle 50%
end
if exist('tempBPM','var')==1 && tempBPM>=min_HR && tempBPM<=max_HR
    startBPM = tempBPM; %replace default with estimated BPM
end

%% HR processing
 
%Isolate max value, calculate bpm and filter anomalies
%For POX, preprocessing normalizes all 'HIGH' values, so max is always the
%1st sample over threshold

ECG_M = zeros(1,N); %preallocate memory
rawHR = [0,0]; %Vector for HR with anomalies
R_counter = 0; %Counter for Nº of beat replacements in signal
T_counter = 0; %Counter for resetting the previous beat
d_change_ratio = change_ratio; %Save default change ratio
HR = [0, startBPM]; %Vector for HR

i = 1; %Global index
while i <= N
    if (N - i) < floor((60/max_HR)*SR)
        L = N-i;
        %Solves infinite loop bug when last value is above threshold
        if L == 0, break; end
    else
        L = floor((60/max_HR)*SR); %Minimum distance between beats)
    end
    %Find the maximum value inside the min_HR window
    if ECG_S(i) > 0
        [m, j] = max(ECG_S(i:i+L));
        %vector to hold only max values
        ECG_M(j+i-1) = m;
        if exist('index1','var') == 0 %Calculate HR only if there's a prev. beat
            index1 = j+i-1; %index of first beat
        else
            index2 = j+i-1; %index of second beat
            %Calculate HR in bpm
            newHR = 60*SR/(index2-index1);
            time = (index2-1)*T;
            %Check that new HR is between min-max range and change ratio
            if newHR <= (HR(end,2)*(1-change_ratio)) ||...
                    newHR >= (HR(end,2)*(1+change_ratio)) ||...
                    newHR < min_HR || newHR > max_HR
                if T_counter >= T_max
                    HR(end+1,1:2) = [time,newHR]; %Reset previous heart beat
                    T_counter = 0;
                else
                    HR(end+1,1:2) = [time,HR(end,2)]; %write previous HR
                    T_counter = T_counter +1;
                    %Verbose replacement time
%                     fprintf('Time:%.1f newHR:%.1f oldHR:%.1f counter:%d\n',...
%                         time, newHR, HR(end,2),T_counter)
                end
                R_counter = R_counter+1;
                change_ratio = d_change_ratio*1.5; %Increase in 50% change ratio
            else
                if T_counter == 1
                    %In case of a skipped beat, average newHR with oldHR
                    HR(end+1,1:2) = [time,(HR(end,2)+newHR)*.5];
                else
                    HR(end+1,1:2) = [time,newHR]; %write new HR
                end
                T_counter = 0;
                change_ratio = d_change_ratio; %change_ratio back to default
            end
            rawHR(end+1,1:2) = [time,newHR];  %write all HR for comparison          
            index1 = index2; %beat 2 now is new beat 1
        end
        i = i+L; %jump minHR distance
    else
        i = i+1;
    end
end         

%Calculate mean HR
meanHR = mean(HR(1:end,2));

%Calculate and print Confidence (ratio of replaced beats vs real beats)
Q = 100-(R_counter/length(HR))*100;
%disp(length(HR))
if length(HR) <= 2
    Q = 0;
    warning('HRcalc:noBeats', 'No beats in ECG signal (check parameters)')
end
if (n_debug==1);fprintf('meanHR %.1f - Confidence %.1f%% \n',meanHR,Q);end

%% Plot
if n_debug == 1
    
    x_limits = [0 (N-1)*T];
    
    %Create an 'error' vector
    errorHR = setdiff(HR,rawHR,'rows');
    
    figure('name','HR Extraction Processing Sequence')
    ax(1) = subplot(4,1,1);
    plot(t,ECG)
    title('Original Signal')
    xlim(x_limits)
    ax(2) = subplot(4,1,2);
    plot(t,ECG_A)
    hold on
    plot(t,S,'r')
    title('Detrend and Absolute values with dynamic threshold value')
    xlim(x_limits)
    ax(3) = subplot(4,1,3);
    stem(t,ECG_S)
    hold on
    stem(t,ECG_M,'r','LineStyle', 'none', 'Marker', '.')
    xlim(x_limits)
    title(sprintf('Values over threshold (mean+R*std) and Peaks, R=%d',R))
    ax(4) = subplot(4,1,4);
    plot(HR(1:end,1),(HR(1:end,2)),'r')
    hold on
    plot(errorHR(1:end,1),(errorHR(1:end,2)),'r','LineStyle', 'none', 'Marker', '.')
    xlim(x_limits)
    ylim([min_HR max_HR])
    title(sprintf('HR - Confidence %.1f%% - * Beat Replacement',Q))
    xlabel('time (s)')
    linkaxes(ax,'x')
end