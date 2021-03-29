clear all
%% Description of task
%
% Your task is to desig your new rate control algorithm for IEEE 802.11.
% For this please use the template 'TransmitRateControlTemplateET4394.m'
% which is based on 'TransmitRateControlExample.m' script from MATLAB WLAN
% Toolbox. Read 'TransmitRateControlTemplateET4394.m' in detail - there you
% will find which parts of the code you should change and where you should
% provide your own algorithm.

% For the evaluation of each project I will look into novelty and
% complexity of the idea. Any constant values that you will use in your
% algorithm should be justified. I will NOT accept the ideas which are
% derivative of the original algorithm provided in 'TransmitRateControlTemplateET4394.m'
% For example changing threshold = [11 14 19 20 25 28 30 31 35] to some
% other vector without specifying WHY and HOW did you change these values
% will not be accepted.

% To gamify the assigment I will evaluate your contribution based on three
% sets of parameters (which I will keep secret and release after deadline). 
% For each set I will randomy change the following parameters ONLY:

% tgacChannel.DelayProfile
% tgacChannel.TransmitReceiveDistance
% meanSNR 
% amplitude
% maxJump
%
% [all these values are explained below in this script]

% The best group will be the one that will have the highest average of:
% (Overall_data_rate_for_scenario_1 + Overall_data_rate_for_scenario_2 +
% Overall_data_rate_for_scenario_3)/3
%
% where 'Overall_data_rate_for_scenario_x' is calculated in line 82 of 
% 'TransmitRateControlTemplateET4394.m',
%
% and the lowest average of:
% (Overall_packet_error_rate_for_scenario_1 +
% Overall_packet_error_rate_for_scenario_2 + 
% Overall_packet_error_rate_for_scenario_3)/3
%
% where 'Overall_packet_error_rate_for_scenario_x' is calculated in line 84 of 
% 'TransmitRateControlTemplateET4394.m',

% This means that to win the game you must design an algorithm that works
% best for all possible contritions.

% Good luck!

% Przemysław Pawełczak, April 1, 2020, 18:15PM


%% Waveform Configuration
%
% Parameters specified by wlanVHTConfig function.
% From here - DO NOT CHANGE THESE VALUES in your simulation

cfgVHT = wlanVHTConfig;   
cfgVHT.ChannelBandwidth = 'CBW40'; % 40 MHz channel bandwidth
cfgVHT.NumUsers = 1;               % Number of users
cfgVHT.UserPositions = 1;          % User position
cfgVHT.NumTransmitAntennas = 1;    % Number of transmit antennas
cfgVHT.PreVHTCyclicShifts = -75;   % Cyclic shift values of additional transmit antennas (not used)
cfgVHT.NumSpaceTimeStreams = 1;    % Number of space-time streams
cfgVHT.SpatialMapping = 'Direct';  % Spatial mapping scheme
cfgVHT.SpatialMappingMatrix = 1;   % Spatial mapping matrix
cfgVHT.Beamforming = 1;            % Enabled beamforming
cfgVHT.STBC = 0;                   % Disabled space time block coding
cfgVHT.MCS = 1;                    % Initial MSC: QPSK rate-1/2
cfgVHT.ChannelCoding = 'BCC';      % FEC coding type is BCC
cfgVHT.APEPLength = 4096;          % APEP length in bytes
cfgVHT.PSDULength;                 % PSDU length (read-only value) 
cfgVHT.GuardInterval = 'Long';     % Guard interval (cyclic prefix) duration
cfgVHT.GroupID = 63;               % Group identification number (default)
cfgVHT.PartialAID = 275;           % Abbreviated indication of PSDU recipients (default)

% To here - DO NOT CHANGE THESE VALUES in your simulation

%% Channel Configuration

tgacChannel = wlanTGacChannel;
sr = wlanSampleRate(cfgVHT);                            % DO NOT CHANGE THIS % Set the sampling rate for the channel
tgacChannel.SampleRate = sr;                            % DO NOT CHANGE THIS % Sample rate (default)                         
tgacChannel.ChannelBandwidth = cfgVHT.ChannelBandwidth; % DO NOT CHANGE THIS % Channel bandwidth
tgacChannel.CarrierFrequency = 5.25e9;                  % DO NOT CHANGE THIS % RF carrier frequency 
tgacChannel.EnvironmentalSpeed = 0.089;                 % DO NOT CHANGE THIS % Speed of the scatterers (default)
tgacChannel.NormalizePathGains = 1;                     % DO NOT CHANGE THIS % Normalize path gains
tgacChannel.UserIndex = 0;                              % DO NOT CHANGE THIS % User index for single or multi-user scenario
tgacChannel.TransmissionDirection = 'Downlink';         % DO NOT CHANGE THIS % Transmission direction
tgacChannel.NumTransmitAntennas = 1;                    % DO NOT CHANGE THIS % Number of transmit antennas
%tgacChannel.TransmitAntennaSpacing = 0.5;               % DO NOT CHANGE THIS % Distance between transmit antenna elements (wavelenght) - only for multi-antenna
tgacChannel.NumReceiveAntennas = 1;                     % DO NOT CHANGE THIS % Number of receive antennas
%tgacChannel.ReceiveAntennaSpacing = 0.5;                % DO NOT CHANGE THIS % Distance between receive antenna elements (wavelenght) - only for multi-antenna
tgacChannel.LargeScaleFadingEffect = 'None';            % DO NOT CHANGE THIS % Large-scale fading effects
%tgacChannel.FluorescentEffect = 1;                      % DO NOT CHANGE THIS % Fluorescent effect
%tgacChannel.PowerLineFrequency = '60Hz';                % DO NOT CHANGE THIS % Power line frequency
tgacChannel.NormalizeChannelOutputs = 1;                % DO NOT CHANGE THIS % Normalize channel outputs
tgacChannel.ChannelFiltering = 1;                       % DO NOT CHANGE THIS % Enable channel filtering
%tgacChannel.NumSamples = 320;                           % DO NOT CHANGE THIS % Used only when ChannelFiltering = 0
%tgacChannel.OutputDataType = 'double';                  % DO NOT CHANGE THIS % Data type of impaired signal
tgacChannel.RandomStream = 'mt19937ar with seed';       % DO NOT CHANGE THIS % Source of random number stream
tgacChannel.Seed = 0;                                   % DO NOT CHANGE THIS % Initial seed of mt19937ar random number stream
tgacChannel.PathGainsOutputPort = 0;                    % DO NOT CHANGE THIS % Enable path gain output

%% Simulation Parameters

numPackets = 1000; % NOTE: the more, the better for relability; Number of packets transmitted during the simulation 

tgacChannel.DelayProfile = 'Model-C';                   % I will select this value myself in final evaluation ---> Delay profile model
tgacChannel.TransmitReceiveDistance = 1;               % I will select this value myself in final evaluation ---> Distance in meters for NLOS
% Select SNR for the simulation
meanSNR = 10;                                                   % I will select this value myself in final evaluation ---> Mean SNR
amplitude = 7;                                                 % I will select this value myself in final evaluation ---> Variation in SNR around the average mean SNR value
baseSNR = sin(linspace(1,10,numPackets)) * amplitude + meanSNR; % DO NOT CHANGE THIS % Generate varying SNR values for each transmitted packet                                                   
maxJump = 3;                                                  % I will select this value myself in final evaluation ---> The maxJump controls the maximum SNR difference between one packet and the next 
snrWalk = baseSNR(1);                                           % DO NOT CHANGE THIS % Set the initial SNR value

%TransmitRateControlTemplateET4394
classifier
% [PUT YOUR OWN RATE CONTROL SCRIP HERE FOLLOWING TransmitRateControlTemplateET4394 TEMPLATE]