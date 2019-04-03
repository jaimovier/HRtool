# HRtool
Matlab tool for processing Heart Rate (HR) signals

#MATLAB usage
Output from 'help HRtool' in command window

HR tool for ECG and POX
  by Javier Jaimovich (2012)
  version 3.4 (updated 21/06/2012)
 
  [HR, meanHR, Q] = HRtool(ECG, SR, debug, par, input)
 
  Returns a matrix with HR (time, bpm) and average HR
  Q is the level of confidence of the extraction (%)
 
  ECG: Input signal (ECG or POX)
  SR: Sampling Rate
  debug: 0 for no debugging, 1 to debug (includes plot)
  par is a vector with [minHR maxHR change_ratio R]
  being:
    minHR: minimum HR limit
    maxHR: maximum HR limit
    R: Ratio between beat and std of ECG for beat threshold (only for ECG)
    change_ratio: Maximum change ratio between beats (e.g. 40% -> 0.4)
    
    Default values:
        debug = 0
        minHR = 50
        maxHR = 140
        change_ratio = 0.3 (30%)
        R = 3
        input = 'ECG'
 
    note: For POX signals HRtool will use 1st sample above threshold for
    calculations. With ECG signals, HRtool will use max absolute value of
    each waveform beat.
    
    ***ISSUES***
    - Q is not considering signals that stop before the end of the file
    (e.g. healthy pox signal for the first 30 sec of a 90 sec file). The
    calculation should be: find de time of the last pulse and if it's more
    than minHR distance from the end of the file, deduct points to Q
    proportionally to that distance (e.g. minHR = 60, last heart beat at 5
    seconds from the end of a 90 sec file, then Q = Q - 100*(5-1)/30; Look
    at eimCorrectHR_Q.m file
    
# Reference
If you are using this tool for academic purposes, please cite the following paper:

Jaimovich, Javier, and R. Benjamin Knapp. 2015. “Creating Biosignal Algorithms for Musical Applications from an Extensive Physiological Database.” In Proceedings of the 2015 Conference on New Interfaces for Musical Expression (NIME 2015). Baton Rouge, LA.

