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
P.SNRRange          = -28:-10;  % SNR Range to simulate in dB
P.NumberTxAntennas  = 2;        % Number of transmission antennas for MIMO
P.NumberRxAntennas  = 3;        % Number of receive antennas for MIMO

% Parameters for AWGN or Bypass channels
% P.ChannelType       = 'AWGN';	% Set 'Bypass' for no channel effect
% P.ChannelLength     = 1;        % It must be one, otherwise error

% Parameters for Multipath channel
P.ChannelType       = 'Multipath';

% P.CDMAUsers = 1;
% P.ChannelLength = 3;
% P.RakeFingers = 3;
% BER1 = simulator(P);
% sim1 = sprintf('Ch. length: %d - Users: %d - Fingers: %d' , P.ChannelLength,P.CDMAUsers,P.RakeFingers);
% 
% P.CDMAUsers = 2;
% P.ChannelLength = 3;
% P.RakeFingers = 3;
% BER2 = simulator(P);
% sim2 = sprintf('Ch. length: %d - Users: %d - Fingers: %d' , P.ChannelLength,P.CDMAUsers,P.RakeFingers);
% 
% P.CDMAUsers = 3;
% P.ChannelLength = 3;
% P.RakeFingers = 3;
% BER3 = simulator(P);
% sim3 = sprintf('Ch. length: %d - Users: %d - Fingers: %d' , P.ChannelLength,P.CDMAUsers,P.RakeFingers);

P.CDMAUsers = 4;
P.ChannelLength = 5;
P.RakeFingers = 5;
BER3 = simulator(P);
sim3 = sprintf('Ch. length: %d - Users: %d - Fingers: %d' , P.ChannelLength,P.CDMAUsers,P.RakeFingers);

figure();
semilogy(P.SNRRange,BER1,'DisplayName',sim1);
hold on;
semilogy(P.SNRRange,BER2,'DisplayName',sim2);
semilogy(P.SNRRange,BER3,'DisplayName',sim3);
xlabel('SNR [dB]','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
xlim([min(P.SNRRange) max(P.SNRRange)]);
grid minor;
legend('-DynamicLegend');