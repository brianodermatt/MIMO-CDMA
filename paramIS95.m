% EPFL Advanced Wireless Receivers
% Project IS95, Spring 2020
% Francesco Gallo, Brian Odermatt

clc; clear all; close all;

% Parameters
P.NumberOfFrames	= 100;      % Total number of sent frames
P.BitsPerUser       = 500;      % Bits per frame that each user is given
P.Modulation        = 1;        % Only BPSK supported!
P.CDMAUsers         = 2;        % Total number of users
P.ConvRate          = 1/2;      % Rate of convolutional code, only 1/2 supported!
P.HamLen            = 64;       % Length of Hadamard Sequence, given in IS95 standard
P.ReceiverType      = 'Rake';	% Only 'Rake' supported!

% Parameters for AWGN or Bypass channels
P.ChannelType       = 'AWGN';	% Set 'Bypass' for no channel effect
P.ChannelLength     = 1;        % It must be one, otherwise error

% % Parameters for Multipath channel
% P.ChannelType       = 'Multipath';
% P.ChannelLength     = 3;

P.SNRRange = -10:20; % SNR Range to simulate in dB

BER = simulator(P);
% BER = simulator_NoDec(P);

simlab = sprintf('%s - Length: %d - Users: %d' ,P.ChannelType,P.ChannelLength,P.CDMAUsers);

figure(1)
semilogy(P.SNRRange,BER,'ro-','DisplayName',simlab,'LineWidth',2)

xlabel('SNR','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
xlim([min(P.SNRRange) max(P.SNRRange)]);
grid minor;
legend('-DynamicLegend');