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
    
    % constellation points for MIMO decoder
    Constellations       = [ 1+1i, 1-1i, -1+1i, -1-1i ];

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
            case 'Bypass'
                H = ones(Users, P.NumberRxAntennas*P.ChannelLength, P.NumberTxAntennas);
            case 'AWGN'
                H = ones(Users, P.NumberRxAntennas*P.ChannelLength, P.NumberTxAntennas);
            case 'Multipath'
                H = sqrt(1/2) * (...
                      randn(Users, P.NumberRxAntennas, P.NumberTxAntennas, P.ChannelLength) +...
                      1i*randn(Users, P.NumberRxAntennas, P.NumberTxAntennas, P.ChannelLength)...
                    );
            otherwise
                error('Channel not supported')
        end

        % Noise initialization (Power = 1 [W])
        snoise = randn(P.NumberRxAntennas,NumOfRXChipsPerUser,Users) + 1i * randn(P.NumberRxAntennas,NumOfRXChipsPerUser,Users);

        % SNR Range
        for ss = 1:length(P.SNRRange)
            SNRdb  = P.SNRRange(ss);
            SNRlin = 10^(SNRdb/10);
            % Normalize noise according to SNR (noise power) and spreading
            % factor (noise is equally distributed over the chips)
            noise  = 1/sqrt(2*SNRlin*SeqLen) * snoise;

            % Channel
            switch P.ChannelType
                case 'Bypass'
                    y = mwaveform;
                case 'AWGN'
                    y = mwaveform + noise;
                case 'Multipath'     
                    y = zeros(P.NumberRxAntennas, NumOfChipsPerUser + P.ChannelLength-1, Users);
                    for i = 1:Users
                       for j = 1:P.NumberRxAntennas
                          for k = 1:P.NumberTxAntennas
                              y(j,:,i) = y(j,:,i) + conv(mwaveform(k,:,i), squeeze(H(i,j,k,:))) + noise(j,:,i);
                          end
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
                rakeAntennas = zeros(P.NumberRxAntennas,P.ChannelLength,NumOfEncBits/Users);                
                for k = 1:P.NumberRxAntennas
                    for l = 1:P.RakeFingers
                        data = y(k,l:l+NumOfChipsPerUser-1,i);
                        rxvecs  = reshape(data,SeqLen,NumOfEncBits/Users);
                        rakeAntennas(k,l,:) = 1/SeqLen * UserSequence.' * rxvecs;
                    end
                end
            
                % Then do MIMO for each of the virtual RAKE antennas before
                % combining them
                sTilde = zeros(P.RakeFingers, P.NumberTxAntennas, NumOfEncBits/Users);
                sHat = zeros(size(sTilde));
                for l = 1:P.RakeFingers
                   Heff = squeeze(H(i,:,:,l));
                   
                   switch P.MIMODetectorType
                        case 'ZeroForcing'
                            H_H = Heff';
                            G = inv(H_H * Heff) * H_H;
                            sTilde(l,:,:) = G * squeeze(rakeAntennas(:,l,:));
                            % could map to closest constellation point but
                            % this would destroy information
%                             for k = 1:P.NumberTxAntennas
%                                 for m = 1:(NumOfEncBits/Users)
%                                    [~,closestIndex] = min(sTilde(l,k,i) - Constellations);
%                                    sHat(l,k,m) = Constellations(closestIndex);
%                                 end
%                             end
                        otherwise
                            error('MIMO Detector not supported');
                    end
                end
                
                % from sTilde, which is available for each RAKE finger on
                % each Tx antenna
                % combine the RAKE fingers
                combinedRake = zeros(P.NumberTxAntennas, NumOfEncBits/Users);
                for k = 1:P.NumberTxAntennas
                    fingersToCombine = squeeze(sTilde(:,k,:));
                    
                    % since we inverted the channel effects already in the
                    % MIMO detection, we can just average the contributions
                    averager = ones(1, P.ChannelLength);
                    combinedRake(k,:) = (1/norm(averager)) * conj(averager) * fingersToCombine;
                end
                
                % since all antennas sent the same data, we can average the
                % estimates of each Tx antenna
                mrc = mean(combinedRake,1);
                
                % Decoding the bits: soft Viterbi decoder
                DecodedBits = step(decoder, real(mrc).');
                
                % Eliminating convolution tails
                RxBits(i,:) = DecodedBits(1:P.BitsPerUser);
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