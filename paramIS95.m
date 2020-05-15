% EPFL Advanced Wireless Receivers
% Project IS95, Spring 2020
% Francesco Gallo, Brian Odermatt

clc; clear all; close all;

% Parameters
P.NumberOfFrames      = 100;
P.NumberOfSymbols     = 1000;

P.AccessType = 'CDMA';
P.CDMAUsers     = 2;

P.Modulation    = 2;        % 2: QPSK

P.ChannelType   = 'Multipath'; % 'AWGN'
P.ChannelLength = 3;

P.HamLen = 64; % Length of Hadamard Sequence, given in IS95 standard

P.SNRRange = -10:20; % SNR Range to simulate in dB

P.ReceiverType  = 'Rake';

BER = simulator(P);

simlab = sprintf('%s - Length: %d - Users: %d' ,P.ChannelType,P.ChannelLength,P.CDMAUsers);

figure(1)
semilogy(P.SNRRange,BER,'ro-','DisplayName',simlab,'LineWidth',2)

xlabel('SNR','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
xlim([min(P.SNRRange) max(P.SNRRange)]);
grid minor;
legend('-DynamicLegend');