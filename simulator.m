% EPFL - Advanced Wireless Receivers
% Final Project:
% CDMA IS95 standard, system simulator
% Brian Odermatt, Francesco Gallo
% May 2020

function BER = simulator(P)

    if P.CDMAUsers > P.HadLen
       error('More users than available sequences')
    end
    
    if P.RakeFingers > P.ChannelLength
       warning('More Rx fingers than available channel taps. Setting the number of fingers to the maximum possible value.')
       P.RakeFingers = P.ChannelLength;
    end
    
    % Initial values and constants
    ConstraintLength     = P.ConstrLen;
    Users                = P.CDMAUsers;
    
    HadamardMatrix       = hadamard(P.HadLen)/sqrt(P.HadLen);
    SpreadSequence       = HadamardMatrix;
    SeqLen               = P.HadLen;
    
    % Total number of bits per frame:
    NumOfBits            = P.BitsPerUser * Users;
    
    % Total number of encoded bits per frame: the convolutional encoder
    % adds a tail of bits that does not carry information.
    NumOfEncBits         = (NumOfBits + (ConstraintLength-1)*Users)/P.ConvRate;
    
    % Number of chips per user, per frame:
    NumOfChipsPerUser    = NumOfEncBits * SeqLen / Users;
    
    % generating PN sequence
    PNSequence           = genPN(NumOfChipsPerUser);

    % Convolutional encoder
    ConvolutionalGeneratorPolynoms = [753 561];
    encoder = comm.ConvolutionalEncoder(...
        'TerminationMethod', 'Terminated',...
        'TrellisStructure', poly2trellis(ConstraintLength, ConvolutionalGeneratorPolynoms)...
    );
    
    % Convolutional decoder
    decoder = comm.ViterbiDecoder(...
        'TerminationMethod', 'Terminated',...
        'TrellisStructure', poly2trellis(ConstraintLength, ConvolutionalGeneratorPolynoms)...
    );
    
    % Channel
    switch P.ChannelType
        case 'Multipath'
            % The convolution of the noisy symbols with channel impulse
            % response adds a convolution tail:
            NumOfRXChipsPerUser = NumOfChipsPerUser+P.ChannelLength-1;
        otherwise
            NumOfRXChipsPerUser = NumOfChipsPerUser;
    end
    
    
    Results = zeros(1,length(P.SNRRange));


    for ii = 1:P.NumberOfFrames

        disp(['Simulating: ' sprintf('%.2f', ii/P.NumberOfFrames*100) '%'])
        
        % Information bits (already initialized in a user matrix to avoid
        % reshaping for the encoder)
        Bits = randi([0 1], NumOfBits/Users, Users); % Random Data

        % Convolutional encoding: rate 1/2
        EncBits = zeros(NumOfEncBits/Users, Users);
        for i = 1:Users
            EncBits(:,i) = step(encoder, Bits(:,i));
        end

        % BPSK Modulation
        symbols = -(2*EncBits - 1);

        % Orthogonal spreading
        txsymbols = SpreadSequence(:,1:Users) * symbols.';

        % Applying PN sequence
        waveform = txsymbols(:).*PNSequence;

        % Reshape to add multi-user antenna suppport
        waveform  = reshape(waveform,1,NumOfChipsPerUser);
        mwaveform = repmat(waveform,[P.NumberTxAntennas 1 Users]);

        % Channel
        switch P.ChannelType
            case 'Multipath'
                H = sqrt(1/2) * (...
                    randn(Users, P.ChannelLength * P.NumberRxAntennas, P.NumberTxAntennas) + ...
                    1i * randn(Users, P.ChannelLength * P.NumberRxAntennas, P.NumberTxAntennas) ...
                );
            otherwise
                error('Channel not supported')
        end

        % Noise initialization (Power = 1 [W]). Independent noise
        % realisation on each (virtual) Rx antenna
        snoise = randn(P.ChannelLength * P.NumberRxAntennas, NumOfRXChipsPerUser, Users) + ...
                 1i * randn(P.ChannelLength * P.NumberRxAntennas, NumOfRXChipsPerUser, Users);

        % SNR Range
        for ss = 1:length(P.SNRRange)
            SNRdb  = P.SNRRange(ss);
            SNRlin = 10^(SNRdb/10);
            
            % Normalize noise according to SNR (noise power) and spreading
            % factor (noise is equally distributed over the chips)
            noise  = 1/sqrt(2*SNRlin*SeqLen) * snoise;

            % Channel
            switch P.ChannelType
                case 'Multipath'
                    y = zeros(P.NumberRxAntennas, NumOfChipsPerUser + P.ChannelLength-1, Users);
                    for i = 1:Users
                        % Reshape into the form of eq. 2 of the assignment
                        HConv = reshape(H(i,:,:), [P.ChannelLength, P.NumberRxAntennas, P.NumberTxAntennas]);
                        for j = 1:P.NumberRxAntennas
                            for k = 1:P.NumberTxAntennas
                                % perform convolution as in eq. 5/6 of the assignment
                                h_jk = HConv(:,j,k);
                                y(j,:,i) = y(j,:,i) + conv(mwaveform(k,:,i), h_jk);
                            end
                            y(j,:,i) = y(j,:,i) + noise(j,:,i);
                        end
                    end
                otherwise
                    error('Channel not supported')
            end
            
            RxBits = zeros(Users,NumOfBits/Users);
            for i = 1:Users
                % first split up into virtual RAKE antennas. There are
                % P.NumberRxAntennas*P.ChannelLength virtual antennas per user
                UserSequence = SpreadSequence(:,i); 
                rakeAntennas = zeros(P.NumberRxAntennas * P.ChannelLength, NumOfEncBits/Users);
                for k = 1:P.NumberRxAntennas
                    for m = 1:P.RakeFingers
                        data = y(k,m:m+NumOfChipsPerUser-1,i) ./ PNSequence.';
                        rxvecs  = reshape(data,SeqLen,NumOfEncBits/Users);
                        rakeAntennas((k-1)*P.RakeFingers + m, :) = UserSequence.' * rxvecs;
                    end
                end
            
                % MIMO with the virtual RAKE antennas directly gives the
                % estimation of the sent signal on each antenna
                H_User = squeeze(H(i,:,:));
                switch P.MIMODetectorType
                    case 'ZF'
                        H_H = H_User';
                        G = inv(H_H * H_User) * H_H;
                        sTilde = G * rakeAntennas;
                        
                    case 'MMSE'
                        H_H = H_User';
                        mult = H_H * H_User;
                        G = inv(mult + P.NumberTxAntennas*SNRlin*eye(size(mult))) * H_H;
                        sTilde = G * rakeAntennas;
                        
                    case 'SIC'
                        yi = rakeAntennas;
                        Hi = H_User;
                        Constellations = [-1 1];    % here, the signal is BPSK
                        for i = 1:P.NumberTxAntennas
                            Hi_H = Hi';
                            Hi_inv = inv(Hi_H * Hi) * Hi_H;
                            g1i_star = Hi_inv(1,:);
                            % here we actually produce sHats and not sTilde
                            for kk = 1:NumOfEncBits/Users
                                temp = g1i_star*yi;
                                [~, closestIndex] = min(temp(kk) - Constellations);
                                sTilde(i,kk) = Constellations(closestIndex);
                            end
                            yi = yi - Hi(:,1) * sTilde(i);
                            Hi(:,1) = [];
                        end
                        
                    otherwise
                        error('MIMO Detector not supported');
                end
                
                % could map to closest constellation point but
                % this would destroy information
                
                % since all antennas sent the same data, we can average the
                % estimates of each Tx antenna
                mrc = mean(sTilde, 1);
                
                % Decoding the bits: soft Viterbi decoder
                decodedBits = step(decoder, real(mrc).');
                
                % Eliminating convolution tails
                RxBits(i,:) = decodedBits(1:P.BitsPerUser);
            end

            % Flatten the bit vectors for BER count
            Bits    = reshape(Bits, NumOfBits, 1);
            RxBits  = reshape(RxBits.', NumOfBits, 1);

            % BER count
            errors      = sum(RxBits ~= Bits);
            Results(ss) = Results(ss) + errors;

        end
    end

    BER = Results/(NumOfBits*P.NumberOfFrames);
end