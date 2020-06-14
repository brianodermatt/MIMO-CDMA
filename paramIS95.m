% EPFL - Advanced Wireless Receivers
% Final Project:
% CDMA IS95 standard, parameter file
% Brian Odermatt, Francesco Gallo
% May 2020

clc; clear all; close all;

% Parameters
P.NumberOfFrames	= 100;      % Total number of sent frames
P.BitsPerUser       = 172;      % Bits per frame that each user is given, for 9600bps frame according to spec 
P.ConvRate          = 1/2;      % Rate of convolutional code, only 1/2
P.ConstrLen         = 9;        % Constraint length of convolutional encoder
P.HadLen            = 64;       % Length of Hadamard Sequence, given in IS95 standard

%P.SNRRange          = -28:-10;  % SNR Range to simulate in dB
SNRRange            = -20:1:0;
P.SNRRange          = SNRRange - 10*log10(P.HadLen);

% P.NumberTxAntennas  = 2;        % Number of transmission antennas for MIMO
% P.NumberRxAntennas  = 2;        % Number of receive antennas for MIMO

% Parameters for Multipath channel
% P.ChannelType       = 'Multipath';   % Since MIMO is implemented, only multipath is possible (inverting an all-1 matrix gives a singularity)

% Parameter for MIMO detection
% P.MIMODetectorType  = 'MMSE';

displaySnrRange = P.SNRRange + 10*log10(P.HadLen);

% P.CDMAUsers = 2;
% P.ChannelLength = 3;
% P.RakeFingers = 3;
% BER = simulator(P);
% sim3 = sprintf('Ch. len.:%d; Users:%d; Fingers:%d',P.ChannelLength,P.CDMAUsers,P.RakeFingers);

% figure();
% semilogy(SNRRange,BER,'DisplayName',sim3);
% xlabel('SNR [dB]','FontSize',12,'FontWeight','bold');
% ylabel('BER','FontSize',12,'FontWeight','bold');
% xlim([min(SNRRange) max(SNRRange)]);
% grid minor;
% legend('-DynamicLegend');


%% first simulation: vary number of users
P.NumberTxAntennas  = 2;        % Number of transmission antennas for MIMO
P.NumberRxAntennas  = 2;        % Number of receive antennas for MIMO
P.ChannelType       = 'Multipath';   % Since MIMO is implemented, only multipath is possible (inverting an all-1 matrix gives a singularity)
P.MIMODetectorType  = 'MMSE';

P.ChannelLength = 3;
P.RakeFingers = 3;
figure();
users = [1 2 3 10 64];
for i = 1:length(users)
   u = users(i);
   P.CDMAUsers = u;
   name = sprintf('Ch.Len.:%d; Users:%d; Fingers:%d; n_T:%d; n_R:%d Detector:%s',P.ChannelLength,P.CDMAUsers,P.RakeFingers,P.NumberTxAntennas,P.NumberRxAntennas,P.MIMODetectorType);
   semilogy(SNRRange,simulator(P),'DisplayName',name);
   hold on;
end
xlabel('SNR [dB]','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
xlim([min(SNRRange) max(SNRRange)]);
grid minor;
legend('-DynamicLegend');