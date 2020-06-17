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

SNRRange            = -20:1:0;  % SNR Range to simulate in dB
P.SNRRange          = SNRRange - 10*log10(P.HadLen);

displaySnrRange = P.SNRRange + 10*log10(P.HadLen);


% %% first simulation: vary number of users
% P.NumberTxAntennas  = 2;        % Number of transmission antennas for MIMO
% P.NumberRxAntennas  = 2;        % Number of receive antennas for MIMO
% P.ChannelType       = 'Multipath';   % Since MIMO is implemented, only multipath is possible (inverting an all-1 matrix gives a singularity)
% P.MIMODetectorType  = 'MMSE';
% P.ChannelLength     = 3;
% P.RakeFingers       = 3;
% figure();
% users = [1 2 3 10 64];
% for i = 1:length(users)
%    u = users(i);
%    P.CDMAUsers = u;
%    name = sprintf('Ch.Len.:%d; Users:%d; Fingers:%d; n_T:%d; n_R:%d Detector:%s',P.ChannelLength,P.CDMAUsers,P.RakeFingers,P.NumberTxAntennas,P.NumberRxAntennas,P.MIMODetectorType);
%    semilogy(SNRRange,simulator(P),'DisplayName',name);
%    hold on;
% end
% xlabel('SNR [dB]','FontSize',12,'FontWeight','bold');
% ylabel('BER','FontSize',12,'FontWeight','bold');
% xlim([min(SNRRange) max(SNRRange)]);
% grid minor;
% legend('-DynamicLegend');
% 
% %% second simulation: vary number of antennas
% P.ChannelType       = 'Multipath';   % Since MIMO is implemented, only multipath is possible (inverting an all-1 matrix gives a singularity)
% P.MIMODetectorType  = 'MMSE';
% P.CDMAUsers         = 2;
% P.ChannelLength     = 3;
% P.RakeFingers       = 3;
% figure();
% antennas = [1,1; 2,2; 4,4; 2,4];  % [n_T, n_R], n_R must be bigger than n_T
% for i = 1:size(antennas,1)
%     nT = antennas(i,1);
%     nR = antennas(i,2);
%     P.NumberTxAntennas  = nT;        % Number of transmission antennas for MIMO
%     P.NumberRxAntennas  = nR;        % Number of receive antennas for MIMO
%     name = sprintf('Ch.Len.:%d; Users:%d; Fingers:%d; n_T:%d; n_R:%d Detector:%s',P.ChannelLength,P.CDMAUsers,P.RakeFingers,P.NumberTxAntennas,P.NumberRxAntennas,P.MIMODetectorType);
%     semilogy(SNRRange,simulator(P),'DisplayName',name);
%     hold on;
% end
% xlabel('SNR [dB]','FontSize',12,'FontWeight','bold');
% ylabel('BER','FontSize',12,'FontWeight','bold');
% xlim([min(SNRRange) max(SNRRange)]);
% grid minor;
% legend('-DynamicLegend');
% 
% %% third simulation: vary channel length
% P.NumberTxAntennas  = 2;        % Number of transmission antennas for MIMO
% P.NumberRxAntennas  = 2;        % Number of receive antennas for MIMO
% P.ChannelType       = 'Multipath';   % Since MIMO is implemented, only multipath is possible (inverting an all-1 matrix gives a singularity)
% P.MIMODetectorType  = 'MMSE';
% P.CDMAUsers         = 2;
% figure();
% lengths = [1, 3, 6];
% for i = 1:length(lengths)
%     channelLength = lengths(i);
%     P.ChannelLength = channelLength;
%     P.RakeFingers = channelLength;
%     name = sprintf('Ch.Len.:%d; Users:%d; Fingers:%d; n_T:%d; n_R:%d Detector:%s',P.ChannelLength,P.CDMAUsers,P.RakeFingers,P.NumberTxAntennas,P.NumberRxAntennas,P.MIMODetectorType);
%     semilogy(SNRRange,simulator(P),'DisplayName',name);
%     hold on;
% end
% xlabel('SNR [dB]','FontSize',12,'FontWeight','bold');
% ylabel('BER','FontSize',12,'FontWeight','bold');
% xlim([min(SNRRange) max(SNRRange)]);
% grid minor;
% legend('-DynamicLegend');
% 
% %% fourth simulation: vary number of RAKE fingers
% P.NumberTxAntennas  = 2;        % Number of transmission antennas for MIMO
% P.NumberRxAntennas  = 2;        % Number of receive antennas for MIMO
% P.ChannelType       = 'Multipath';   % Since MIMO is implemented, only multipath is possible (inverting an all-1 matrix gives a singularity)
% P.MIMODetectorType  = 'MMSE';
% P.ChannelLength     = 4;
% P.CDMAUsers         = 2;
% figure();
% fingers = [1, 2, 3, 4];
% for i = 1:length(fingers)
%     finger = fingers(i);
%     P.RakeFingers = finger;
%     name = sprintf('Ch.Len.:%d; Users:%d; Fingers:%d; n_T:%d; n_R:%d Detector:%s',P.ChannelLength,P.CDMAUsers,P.RakeFingers,P.NumberTxAntennas,P.NumberRxAntennas,P.MIMODetectorType);
%     semilogy(SNRRange,simulator(P),'DisplayName',name);
%     hold on;
% end
% xlabel('SNR [dB]','FontSize',12,'FontWeight','bold');
% ylabel('BER','FontSize',12,'FontWeight','bold');
% xlim([min(SNRRange) max(SNRRange)]);
% grid minor;
% legend('-DynamicLegend');

%% fifth simulation: different detectors
P.NumberTxAntennas  = 2;        % Number of transmission antennas for MIMO
P.NumberRxAntennas  = 2;        % Number of receive antennas for MIMO
P.ChannelType       = 'Multipath';   % Since MIMO is implemented, only multipath is possible (inverting an all-1 matrix gives a singularity)
P.ChannelLength     = 3;
P.RakeFingers       = 3;
P.CDMAUsers         = 2;
figure();
detectors = ["ZF", "MMSE", "SIC"];
%detectors = ["SIC"];

for i = 1:length(detectors)
    detector = detectors(i);
    P.MIMODetectorType = detector;
    name = sprintf('Ch.Len.:%d; Users:%d; Fingers:%d; n_T:%d; n_R:%d Detector:%s',P.ChannelLength,P.CDMAUsers,P.RakeFingers,P.NumberTxAntennas,P.NumberRxAntennas,P.MIMODetectorType);
    semilogy(SNRRange,simulator(P),'DisplayName',name);
    hold on;
end
xlabel('SNR [dB]','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
xlim([min(SNRRange) max(SNRRange)]);
grid minor;
legend('-DynamicLegend');