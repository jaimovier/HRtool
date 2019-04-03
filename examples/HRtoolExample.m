%% HRtool example
% by Javier Jaimovich

clc; help HRtool.m

%Load examples with ECG or POX signals
load('HRtoolExamples.mat')

%Fill array with names
for i=1:size(Examples,1)
    MenuItems(i,1) = {eval(sprintf('Examples{%d}.name',i))};
end

%Assign selected ECG signal to ECG var
choice = menu('Select Example Signal',MenuItems);

signal = Examples{choice}.signal;
SR = Examples{choice}.SR;
type = Examples{choice}.type;

%print signal's details
disp(Examples{choice})

%% Run HRtool

%Set HRtool parameters
minHR = 50;
maxHR = 150;
ChangeRatio = 0.2;
Threshold = 2;

par = [minHR maxHR ChangeRatio Threshold]; %make par vector
[HR, meanHR, Q] = HRtool(signal,SR,1,par,type);