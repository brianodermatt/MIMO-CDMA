% EPFL - Advanced Wireless Receivers
% Final Project:
% CDMA IS95 standard, parameter file
% Brian Odermatt, Francesco Gallo
% May 2020

clc; clear; close all;

% Fixed Parameters
P.NumberOfFrames	= 100;      % Total number of sent frames
P.BitsPerUser       = 172;      % Bits per frame that each user is given, for 9600bps frame according to spec 
P.HadLen            = 64;       % Length of Hadamard Sequence, given in IS95 standard
P.ConvRate          = 1/2;      % Rate of convolutional code, only 1/2
P.ConstrLen         = 9;        % Constraint length of convolutional encoder
P.SNRRange          = -10:1:10;  % SNR Range to simulate in dB
SNR                 = P.SNRRange;

%% Parameters for Bypass or AWGN Simulation
% P.NumberTxAntennas  = 3;        % Number of transmitter antennas
% P.NumberRxAntennas  = 3;        % Number of receiver antennas
% P.ChannelType       = 'AWGN';
% P.ChannelLength     = 1;
% P.RakeFingers       = 1;
% P.CDMAUsers = 2;
% BER = simulator(P);

%% first simulation: vary number of users
P.NumberTxAntennas  = 2;        % Number of transmitter antennas
P.NumberRxAntennas  = 2;        % Number of receiver antennas
P.ChannelType       = 'Multipath';
P.MIMODetectorType  = 'MMSE';
P.ChannelLength     = 3;
P.RakeFingers       = 3;
users = [1 16 32 64];

figure();
for i = 1:length(users)
   P.CDMAUsers = users(i);
   semilogy(SNR, simulator(P), 'DisplayName', sprintf('%d users', users(i)));
   hold on;
end
title(sprintf('MIMO CDMA: n_T=%d, n_R=%d, %s, L_C=%d, N_R=%d, %s', P.NumberTxAntennas, P.NumberRxAntennas, P.ChannelType, P.ChannelLength, P.RakeFingers, P.MIMODetectorType));
xlabel('SNR [dB]','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
grid minor;
legend('-DynamicLegend', 'Location','southwest');

%% second simulation: vary number of antennas
P.ChannelType       = 'Multipath';
P.MIMODetectorType  = 'MMSE';
P.ChannelLength     = 3;
P.RakeFingers       = 3;
P.CDMAUsers         = 2;
antennas = [1,1; 2,2; 4,4; 2,4];  % [N_TX, N_RX], N_RX >= N_TX

figure();
for i = 1:size(antennas,1)
    P.NumberTxAntennas = antennas(i,1);    % Number of transmitter antennas
    P.NumberRxAntennas = antennas(i,2);    % Number of receiver antennas
    semilogy(SNR, simulator(P), 'DisplayName', sprintf('n_T = %d; n_R = %d', antennas(i,1), antennas(i,2)));
    hold on;
end
title(sprintf('MIMO CDMA: %s, L_C=%d, N_R=%d, %s, %d users', P.ChannelType, P.ChannelLength, P.RakeFingers, P.MIMODetectorType, P.CDMAUsers));
xlabel('SNR [dB]','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
grid minor;
legend('-DynamicLegend', 'Location','southwest');

%% third simulation: vary channel length
P.NumberTxAntennas  = 2;        % Number of transmitter antennas
P.NumberRxAntennas  = 2;        % Number of receiver antennas
P.ChannelType       = 'Multipath';
P.MIMODetectorType  = 'MMSE';
P.CDMAUsers         = 2;
taps = [1, 3, 6];

figure();
for i = 1:length(taps)
    P.ChannelLength = taps(i);
    P.RakeFingers   = taps(i);
    semilogy(SNR, simulator(P), 'DisplayName', sprintf('Channel length %d', taps(i)));
    hold on;
end

title(sprintf('MIMO CDMA: n_T=%d, n_R=%d, %s, max. number of RAKE fingers, %s, %d users', P.NumberTxAntennas, P.NumberRxAntennas, P.ChannelType, P.MIMODetectorType, P.CDMAUsers));
xlabel('SNR [dB]','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
grid minor;
legend('-DynamicLegend', 'Location','southwest');

%% fourth simulation: vary number of RAKE fingers
P.NumberTxAntennas  = 2;        % Number of transmitter antennas
P.NumberRxAntennas  = 2;        % Number of receiver antennas
P.ChannelType       = 'Multipath';
P.MIMODetectorType  = 'MMSE';
P.ChannelLength     = 4;
P.CDMAUsers         = 2;
fingers = [1, 2, 3, 4];

figure();
for i = 1:length(fingers)
    P.RakeFingers = fingers(i);
    semilogy(SNR, simulator(P), 'DisplayName', sprintf('%d RAKE fingers', fingers(i)));
    hold on;
end

title(sprintf('MIMO CDMA: n_T=%d, n_R=%d, %s, L_C=%d, %s, %d users', P.NumberTxAntennas, P.NumberRxAntennas, P.ChannelType, P.ChannelLength, P.MIMODetectorType, P.CDMAUsers));
xlabel('SNR [dB]','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
grid minor;
legend('-DynamicLegend', 'Location','southwest');

%% fifth simulation: different detectors
P.NumberTxAntennas  = 2;        % Number of transmitter antennas
P.NumberRxAntennas  = 2;        % Number of receiver antennas
P.ChannelType       = 'Multipath';
P.ChannelLength     = 3;
P.RakeFingers       = 3;
P.CDMAUsers         = 2;

detectors = {'ZF'; 'MMSE'; 'SIC'};

figure();
for i = 1:length(detectors)
    P.MIMODetectorType = detectors{i};
    semilogy(SNR, simulator(P), 'DisplayName', sprintf('%s detector', detectors{i}));
    hold on;
end
title('MIMO CDMA: N_{TX} = N_{RX} = 2, Multipath, 3 taps, 3 rake f., 2 users')
title(sprintf('MIMO CDMA: n_T=%d, n_R=%d, %s, L_C=%d, N_R=%d, %d users', P.NumberTxAntennas, P.NumberRxAntennas, P.ChannelType, P.ChannelLength, P.RakeFingers, P.CDMAUsers));
xlabel('SNR [dB]','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
grid minor;
legend('-DynamicLegend', 'Location','southwest');