% EPFL - Advanced Wireless Receivers
% Final Project:
% CDMA IS95 standard, system simulator
% Brian Odermatt, Francesco Gallo
% May 2020


% Indices in the code:
% To improve readability, any two or more for loops that run over the same
% variable have the same indices 
% ii ----> Runs over the users
% jj ----> Runs over the Rx antennas
% kk ----> Runs over the Tx antennas
% mm ----> Runs over the Rake fingers
% bb ----> Runs over the bits of a string


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
    SpreadSequence       = hadamard(P.HadLen)/sqrt(P.HadLen);
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
        'TerminationMethod', 'Terminated', 'TrellisStructure', ...
         poly2trellis(ConstraintLength, ConvolutionalGeneratorPolynoms) );
    
    % Convolutional decoder
    decoder = comm.ViterbiDecoder(...
        'TerminationMethod', 'Terminated', 'TrellisStructure', ...
         poly2trellis(ConstraintLength, ConvolutionalGeneratorPolynoms) );
    
    % Channel
    switch P.ChannelType
        case 'Multipath'
            % The convolution of the noisy symbols with channel impulse
            % response adds a convolution tail:
            NumOfRxChipsPerUser = NumOfChipsPerUser+P.ChannelLength-1;
        otherwise
            NumOfRxChipsPerUser = NumOfChipsPerUser;
    end
    
    % Initialize variable for error counting
    Results = zeros(1,length(P.SNRRange));

    for n = 1:P.NumberOfFrames

        disp(['Simulating: ' sprintf('%.2f', n/P.NumberOfFrames*100) '%'])
        
        %% Transmitter
        
        % Information bits: random bit strings in a
        % user matrix: each column represents a user string
        Bits = randi([0 1], NumOfBits/Users, Users);

        % Convolutional encoding: rate 1/2
        EncBits = zeros(NumOfEncBits/Users, Users);
        for ii = 1:Users
            EncBits(:,ii) = step(encoder, Bits(:,ii));
        end

        % BPSK Modulation:
        % 0 ----> +1
        % 1 ----> -1
        symbols = 1 - 2*EncBits;

        % Orthogonal spreading
        % dimensions of txsymbols:
        % (SeqLen x Users) x (Users x NumOfEncBits/Users) = 
        %          SeqLen  x  NumOfEncBits/Users =
        txsymbols = SpreadSequence(:, 1:Users) * symbols.';
        
        % Applying PN sequence
        % txsymbols(:) reshapes all elements of A into a single column 
        % vector, acting columnwise. Dimension:
        % (SeqLen*NumOfEncBits/Users) x 1 = NumOfChipsPerUser x 1
        % PNSequence   has dimension:       NumOfChipsPerUser x 1
        % So, waveform has dimension:       NumOfChipsPerUser x 1
        waveform = txsymbols(:).*PNSequence;

       
        % The following operation is equivalent to the (non-conjugate)
        % transpose operation:
        waveform  = reshape(waveform, 1, NumOfChipsPerUser);
        
        % Add multi-user and multi-TxAntenna suppport
        % Dimensions are now:
        % NumberTxAntennas x NumOfChipsPerUser x Users
        mwaveform = repmat(waveform,[P.NumberTxAntennas 1 Users]);

         
        %% Channel and Noise Realizations

        switch P.ChannelType
            case 'Bypass'
                error('Channel not supported yet')

            case 'AWGN'
                error('Channel not supported yet')

            case 'Singlepath Fading'
                error('Channel not supported yet')
            
            case 'Multipath'
                % MIMO multipath channel matrix 
                H = sqrt(1/2) * (...
                    randn(P.ChannelLength * P.NumberRxAntennas, P.NumberTxAntennas, Users) + ...
                    1i * randn(P.ChannelLength * P.NumberRxAntennas, P.NumberTxAntennas, Users));
            otherwise
                error('Channel not supported')
        end

        % Noise initialization: 
        % independent gaussian entries with unit variance (unit average power)
        snoise = randn(P.ChannelLength * P.NumberRxAntennas, NumOfRxChipsPerUser, Users) + ...
                 1i * randn(P.ChannelLength * P.NumberRxAntennas, NumOfRxChipsPerUser, Users);

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
                
                case 'Bypass'
                    error('Channel not supported yet')

                case 'AWGN'
                    error('Channel not supported yet')

                case 'Singlepath Fading'
                    error('Channel not supported yet')
                
                case 'Multipath'
                    % Add here explanation on dimensions
                    y = zeros(P.NumberRxAntennas, NumOfRxChipsPerUser, Users);
                    
                    for ii = 1:Users
                        % Reshape to consider physical antennas
                        HConv = reshape(H(:,:,ii), [P.ChannelLength, P.NumberRxAntennas, P.NumberTxAntennas]);
                        
                        for jj = 1:P.NumberRxAntennas
                            
                            % perform MIMO multipath channel convolution
                            for kk = 1:P.NumberTxAntennas
                                h_jk = HConv(:,jj,kk);
                                y(jj,:,ii) = y(jj,:,ii) + conv(mwaveform(kk,:,ii), h_jk);
                            end
                            
                            % add noise
                            y(jj,:,ii) = y(jj,:,ii) + noise(jj,:,ii);
                            
                        end
                    end
                    
                otherwise
                    error('Channel not supported')
            end
            
            
            %% Receiver
            
            % Initialize received bits in a user matrix: 
            % each row represents a user string
            RxBits = zeros(Users, NumOfBits/Users);
            
            for ii = 1:Users
                % first split up into virtual RAKE antennas. There are
                % P.NumberRxAntennas*P.ChannelLength virtual antennas per user
                UserSequence = SpreadSequence(:,ii) * SeqLen; 
                    % NOTE: Despreading should bring the signal back to
                    % unit power (in a bypass situation).
                    % Here the rx signal power without SeqLen coefficient gets
                    % decreased by a total factor of 64, that's why I
                    % would multiply by 64 to get again a total average
                    % power of 1 for the signal
                rakeAntennas = zeros(P.NumberRxAntennas * P.ChannelLength, NumOfEncBits/Users);
                
                for jj = 1:P.NumberRxAntennas
                    for mm = 1:P.RakeFingers
                        
                        % Remove PN sequence
                        % data has dimensions: NumOfChipsPerUser x 1
                        data = y(jj, mm:mm+NumOfChipsPerUser-1, ii) ./ PNSequence.';
                        
                        % Reshape to apply despreading 
                        rxvecs  = reshape(data, [SeqLen, NumOfEncBits/Users]);
                        
                        % Orthogonal despreading:
                        % each despreading operation gives rise to a 
                        % (1 x NumOfEncBits/Users) vector
                        rakeAntennas((jj-1)*P.RakeFingers + mm, :) = UserSequence.' * rxvecs;
                        
                    end
                end
            
                % MIMO with the virtual RAKE antennas directly gives the
                % estimate of the sent signal on each antenna
                % squeeze removes dimensions of length 1
                H_User = squeeze(H(:,:,ii));
                
                % MIMO detector
                switch P.MIMODetectorType
                    case 'ZF'
                        G = (H_User' * H_User) \ H_User';
                        sTilde = G * rakeAntennas;
                        
                    case 'MMSE'
                        A = H_User' * H_User;
                        B = P.NumberTxAntennas*SNRlin*eye(size(A));
                        G = (A + B) \ H_User';
                        sTilde = G * rakeAntennas;
                        
                    case 'SIC'
                        yi = rakeAntennas;
                        Hi = H_User;
                        % BPSK modulation
                        Constellations = [-1 1];   
                        for kk = 1:P.NumberTxAntennas
                            Gi = (Hi' * Hi) \ Hi';
                            g1i_star = Gi(1,:);
                            % here we actually produce sHats and not sTilde
                            for bb = 1:NumOfEncBits/Users
                                temp = g1i_star*yi;
                                [~, closestIndex] = min(temp(bb) - Constellations);
                                sTilde(kk,bb) = Constellations(closestIndex);
                            end
                            yi = yi - Hi(:,1) * sTilde(kk);
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
                RxBits(ii,:) = decodedBits(1:P.BitsPerUser);
            end

            % Flatten the bit vectors for BER count
            Bits    = reshape(Bits,     NumOfBits, 1);
            RxBits  = reshape(RxBits.', NumOfBits, 1);

            % Error count
            errors      = sum(RxBits ~= Bits);
            % Add to the errors found in the previous frames
            Results(ss) = Results(ss) + errors;

        end
    end
    
    % Calculate the BER
    BER = Results/(NumOfBits*P.NumberOfFrames);
end