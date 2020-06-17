% EPFL - Advanced Wireless Receivers
% Final Project:
% CDMA IS95 standard, parameter file
% Brian Odermatt, Francesco Gallo
% May 2020

clc; clear all; close all;

% Fixed Parameters
P.NumberOfFrames	= 100;      % Total number of sent frames
P.BitsPerUser       = 172;      % Bits per frame that each user is given, for 9600bps frame according to spec 
P.ConvRate          = 1/2;      % Rate of convolutional code, only 1/2
P.ConstrLen         = 9;        % Constraint length of convolutional encoder
P.HadLen            = 64;       % Length of Hadamard Sequence, given in IS95 standard
P.SNRRange          = -43:-8;   % SNR Range to simulate in dB
EbNoRange           = P.SNRRange + 10*log10(P.HadLen);

%% Parameters for Bypass or AWGN Simulation
% P.NumberTxAntennas  = 3;        % Number of transmitter antennas
% P.NumberRxAntennas  = 3;        % Number of receiver antennas
% P.ChannelType       = 'AWGN';
% P.ChannelLength     = 1;
% P.RakeFingers       = 1;
% P.CDMAUsers = 2;
% BER = simulator(P)

%% first simulation: vary number of users
% P.NumberTxAntennas  = 2;        % Number of transmitter antennas
% P.NumberRxAntennas  = 2;        % Number of receiver antennas
% P.ChannelType       = 'Multipath';
% P.MIMODetectorType  = 'MMSE';
% P.ChannelLength     = 3;
% P.RakeFingers       = 3;
% users = [1 2 8 32 64];

% str = {};
% figure()
% for i = 1:length(users)
%    P.CDMAUsers = users(i);
%    str = [str, strcat(num2str(users(i)), 'users')];
%    semilogy(EbNoRange, simulator(P));
%    hold on;
% end
% title('MIMO CDMA: N_{TX} = N_{RX} = 2, Multipath, 3 taps, 3 rake fing., MMSE')
% xlabel('E_B/N_0 [dB]','FontSize',12,'FontWeight','bold');
% ylabel('BER','FontSize',12,'FontWeight','bold');
% grid minor;
% legend(str{:}, 'Location','southwest')

%% second simulation: vary number of antennas
% P.ChannelType       = 'Multipath';
% P.MIMODetectorType  = 'MMSE';
% P.ChannelLength     = 3;
% P.RakeFingers       = 3;
% P.CDMAUsers         = 2;
% antennas = [1,1; 2,2; 4,4; 2,4];  % [N_TX, N_RX], N_RX >= N_TX
% 
% str = {};
% figure()
% for i = 1:size(antennas,1)
%     P.NumberTxAntennas = antennas(i,1);    % Number of transmitter antennas
%     P.NumberRxAntennas = antennas(i,2);    % Number of receiver antennas
%     str = [str, strcat('N_{TX} = ', num2str(antennas(i,1)), '; N_{RX} = ', num2str(antennas(i,2)))];
%     semilogy(EbNoRange, simulator(P));
%     hold on;
% end
% title('MIMO CDMA: Multipath, 3 taps, 3 rake fing., MMSE, 2 users')
% xlabel('E_B/N_0 [dB]','FontSize',12,'FontWeight','bold');
% ylabel('BER','FontSize',12,'FontWeight','bold');
% grid minor;
% legend(str{:}, 'Location','southwest')

%% third simulation: vary channel length
% P.NumberTxAntennas  = 2;        % Number of transmitter antennas
% P.NumberRxAntennas  = 2;        % Number of receiver antennas
% P.ChannelType       = 'Multipath';
% P.MIMODetectorType  = 'MMSE';
% P.CDMAUsers         = 2;
% taps = [1, 3, 6];

% str = {};
% figure()
% for i = 1:length(taps)
%     P.ChannelLength = taps(i);
%     P.RakeFingers   = taps(i);
%     str = [str, strcat(num2str(taps(i)), ' channel taps')];
%     semilogy(EbNoRange, simulator(P));
%     hold on;
% end
% 
% title('MIMO CDMA: N_{TX} = N_{RX} = 2, Multipath, max. rake fing., MMSE, 2 users')
% xlabel('E_B/N_0 [dB]','FontSize',12,'FontWeight','bold');
% ylabel('BER','FontSize',12,'FontWeight','bold');
% grid minor;
% legend(str{:}, 'Location','southwest')

%% fourth simulation: vary number of RAKE fingers
% P.NumberTxAntennas  = 2;        % Number of transmitter antennas
% P.NumberRxAntennas  = 2;        % Number of receiver antennas
% P.ChannelType       = 'Multipath';
% P.MIMODetectorType  = 'MMSE';
% P.ChannelLength     = 4;
% P.CDMAUsers         = 2;
% fingers = [1, 2, 3, 4];
% 
% str = {};
% figure()
% for i = 1:length(fingers)
%     P.RakeFingers = fingers(i);
%     str = [str, strcat(num2str(fingers(i)), ' rake fingers')];
%     semilogy(EbNoRange, simulator(P));
%     hold on;
% end
% 
% title('MIMO CDMA: N_{TX} = N_{RX} = 2, Multipath, 4 taps, MMSE, 2 users')
% xlabel('E_B/N_0 [dB]','FontSize',12,'FontWeight','bold');
% ylabel('BER','FontSize',12,'FontWeight','bold');
% grid minor;
% legend(str{:}, 'Location','southwest')

%% fifth simulation: different detectors
P.NumberTxAntennas  = 2;        % Number of transmitter antennas
P.NumberRxAntennas  = 2;        % Number of receiver antennas
P.ChannelType       = 'Multipath';
P.ChannelLength     = 3;
P.RakeFingers       = 3;
P.CDMAUsers         = 2;
detectors = ["ZF", "MMSE", "SIC"];

str = {};
figure()
for i = 1:length(detectors)
    P.MIMODetectorType = detectors(i);
    str = [str, strcat(detectors(i), ' detector')];
    semilogy(EbNoRange, simulator(P));
    hold on;
end
title('MIMO CDMA: N_{TX} = N_{RX} = 2, Multipath, 3 taps, 3 rake f., 2 users')
xlabel('E_B/N_0 [dB]','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
grid minor;
legend(str{:}, 'Location','southwest')