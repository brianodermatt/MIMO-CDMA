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
       warning('More Rx fingers than available channel taps.')
       print('Setting the number of fingers to the maximum possible value.')
       P.RakeFingers = P.ChannelLength;
    end
    
    % Initial values and constants
    ConstraintLength     = P.ConstrLen;
    Users                = P.CDMAUsers;
    
    % Normalized Hadamard matrix: 
    % normalize so that the total symbol power does not change
    HadamardMatrix       = hadamard(P.HadLen)/sqrt(P.HadLen);
    SpreadSequence       = HadamardMatrix;
    SeqLen               = P.HadLen;
    
    % Total number of bits per frame:
    NumOfBits            = P.BitsPerUser * Users;
    
    % Total number of encoded bits per frame: 
    % the convolutional encoder adds a termination tail of bits that 
    % does not carry information.
    NumOfEncBits         = (NumOfBits + (ConstraintLength-1)*Users)/P.ConvRate;
    
    % Number of chips per user, per frame:
    NumOfChipsPerUser    = NumOfEncBits * SeqLen / Users;
    
    % Generating PN sequence
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
    
    % Initialize variable for error counting
    Results = zeros(1,length(P.SNRRange));

    for ii = 1:P.NumberOfFrames

        disp(['Simulating: ' sprintf('%.2f', ii/P.NumberOfFrames*100) '%'])
        
        %% Transmitter
        
        % Information bits: random bit strings in a
        % user matrix: each column represents a user string
        Bits = randi([0 1], NumOfBits/Users, Users);

        % Convolutional encoding: rate 1/2
        EncBits = zeros(NumOfEncBits/Users, Users);
        for i = 1:Users
            EncBits(:,i) = step(encoder, Bits(:,i));
        end

        % BPSK Modulation:
        % 0 ----> +1
        % 1 ----> -1
        symbols = 1 - 2*EncBits;

        % Orthogonal spreading
        % Add here explanation on dimensions
        txsymbols = SpreadSequence(:, 1:Users) * symbols.';
        
        % Applying PN sequence
        % Add here explanation on dimensions: elementwise multiplication
        waveform = txsymbols(:).*PNSequence;

        % Reshape to add multi-user antenna suppport:
        % Add here explanation on dimensions
        waveform  = reshape(waveform, 1, NumOfChipsPerUser);
        mwaveform = repmat(waveform,[P.NumberTxAntennas 1 Users]);

        
        %% Channel and Noise Realizations

        switch P.ChannelType
            case 'Multipath'
                % MIMO multipath channel matrix 
                % Add here explanation on dimensions
                % Add explanation on distribution of channel coefficients
                H = sqrt(1/2) * (...
                    randn(P.ChannelLength * P.NumberRxAntennas, P.NumberTxAntennas, Users) + ...
                    1i * randn(P.ChannelLength * P.NumberRxAntennas, P.NumberTxAntennas, Users));
            otherwise
                error('Channel not supported')
        end

        % Noise initialization: 
        % independent gaussian entries with unit variance (unit average power)
        % Add here explanation on dimensions
        snoise = randn(P.ChannelLength * P.NumberRxAntennas, NumOfRXChipsPerUser, Users) + ...
                 1i * randn(P.ChannelLength * P.NumberRxAntennas, NumOfRXChipsPerUser, Users);

        % SNR Range
        for ss = 1:length(P.SNRRange)
            SNRdb  = P.SNRRange(ss);
            SNRlin = 10^(SNRdb/10);
            
            % Normalize noise according to 
            % 1. SNR: If signal is normalized to 1, 1/SNR is the noise power;
            % 2. Factor 2: Noise power is equally distributed over the real and imaginary parts;
            % 3. Spread factor: noise is equally distributed over the chips.
            norm_factor = 1/sqrt(2*SNRlin*SeqLen);
            noise  = norm_factor * snoise;

            %% Channel Transmission
            switch P.ChannelType
                case 'Multipath'
                    
                    % Add here explanation on dimensions
                    y = zeros(P.NumberRxAntennas, NumOfChipsPerUser+P.ChannelLength-1, Users);
                    
                    for i = 1:Users
                        % Reshape into the form of eq. 2 of the assignment
                        % Add here explanation on dimensions
                        HConv = reshape(H(:,:,i), [P.ChannelLength, P.NumberRxAntennas, P.NumberTxAntennas]);
                        
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
            
            
            %% Receiver
            
            % Initialize received bits in a user matrix: 
            % each row represents a user string
            RxBits = zeros(Users, NumOfBits/Users);
            
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
                H_User = squeeze(H(:,:,i));
                
                % MIMO detector
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