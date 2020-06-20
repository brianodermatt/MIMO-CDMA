% EPFL - Advanced Wireless Receivers
% Final Project:
% CDMA IS95 standard, system simulator
% Brian Odermatt, Francesco Gallo
% May 2020


% Indices in the code:
% To improve readability, any loop that run over a specific variable has a 
% specific index
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
    
    % Hadamard matrix: 
    % normalization is performed on the noise: noise power is multiplied by 64
    SpreadSequence       = hadamard(P.HadLen);
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
    % PN sequence is not normalized, the noise power will be therefore
    % increased by the following factor
    PNPower              = 2;

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
        
        % Flatten the bit vectors for BER count
        Bits = reshape(Bits, NumOfBits, 1);

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
        % txsymbols(:) reshapes all elements into a single column 
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
                P.ChannelLength = 1;
                P.RakeFingers   = 1;
                P.MIMODetectorType = 'Simple';
                
            case 'AWGN'
                P.ChannelLength = 1;
                P.RakeFingers   = 1;
                P.MIMODetectorType = 'Simple';
            
            case 'Multipath'
                % MIMO multipath channel matrix
                H = 1/sqrt(2) * (randn(P.ChannelLength * P.NumberRxAntennas, P.NumberTxAntennas, Users) + ...
                    1/sqrt(2) * 1i * randn(P.ChannelLength * P.NumberRxAntennas, P.NumberTxAntennas, Users));
            otherwise
                error('Channel not supported')
        end

        % Noise initialization: 
        % independent gaussian entries with unit variance (unit average power)
        % Noise is added to each chip after the convolution of the symbols 
        % with the channel impulse response
        snoise = 1/sqrt(2) * randn(P.ChannelLength * P.NumberRxAntennas, NumOfRxChipsPerUser, Users) + ...
                 1/sqrt(2) * 1i * randn(P.ChannelLength * P.NumberRxAntennas, NumOfRxChipsPerUser, Users);

        % SNR Range
        for ss = 1:length(P.SNRRange)
            SNRdb  = P.SNRRange(ss);
            SNRlin = 10^(SNRdb/10);
            
            % To allow fair simulation, noise needs to be multiplied to
            % compensate for effects of spreading and code rate:
            % 1/SNR is the noise power, if defined for a signal power of 1
            % Factor 1/2: noise power is equally split over the real and imaginary parts
            % ConvRate: Convolutional encoder of rate 1/2. 2 encoded bits
            %   for every information bit.
            % Spread factor (SeqLen): noise is equally distributed over the
            %   chips; since the chips are not normalized, we increase
            %   noise power.
            % PNPower: PN sequence is not normalized, therefore also the
            %   noise power must be increased accordingly
            % Channel length: The more taps, the more power.
            % Tx Antennas: The more antennas, the more power.
            %   Instead of normalizing the channel taps with the channel
            %   length and the number of transmitting antennas, we do it
            %   increasing noise power
            normFactor = sqrt(1/2 * 1/SNRlin * 1/P.ConvRate * SeqLen * PNPower * P.NumberTxAntennas * P.ChannelLength);
            noise = normFactor * snoise;

            %% Channel Transmission
            
            switch P.ChannelType
                
                case 'Bypass'
                    y = zeros(P.NumberRxAntennas, NumOfRxChipsPerUser, Users);
                    for ii = 1:Users
                        for jj = 1:P.NumberRxAntennas
                            for kk = 1:P.NumberTxAntennas
                                y(jj,:,ii) = y(jj,:,ii) + mwaveform(kk,:,ii);
                            end                            
                        end
                    end
                    y = 1/sqrt(P.NumberTxAntennas) * y;

                case 'AWGN'
                    y = zeros(P.NumberRxAntennas, NumOfRxChipsPerUser, Users);
                    for ii = 1:Users
                        for jj = 1:P.NumberRxAntennas
                            for kk = 1:P.NumberTxAntennas
                                y(jj,:,ii) = y(jj,:,ii) + mwaveform(kk,:,ii);
                            end
                            % add noise
                            y(jj,:,ii) = y(jj,:,ii) + noise(jj,:,ii);
                        end
                    end
                
                case 'Multipath'
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
            RxBits = zeros(NumOfBits/Users, Users);
            
            for ii = 1:Users
                % first split up into virtual RAKE antennas. There are
                % P.NumberRxAntennas*P.ChannelLength virtual antennas per user
                UserSequence = SpreadSequence(:,ii); 

                VirtualAntennas = zeros(P.NumberRxAntennas * P.ChannelLength, NumOfEncBits/Users);
                
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
                        VirtualAntennas((jj-1)*P.RakeFingers + mm, :) = UserSequence.' * rxvecs;
                    end
                end
            
                % MIMO with the virtual RAKE antennas directly gives the
                % estimate of the sent signal on each antenna
                % squeeze removes dimensions of length 1
                
                % MIMO detector
                switch P.MIMODetectorType
                    
                    case 'Simple'
                        % For AWGN and Bypass channels
                        sTilde = VirtualAntennas;
                    
                    case 'ZF'
                        H_ii = squeeze(H(:,:,ii));
                        G = (H_ii' * H_ii) \ H_ii';
                        sTilde = G * VirtualAntennas;
                        
                    case 'MMSE'
                        H_ii = squeeze(H(:,:,ii));
                        A = H_ii' * H_ii;
                        % noise variance is given by 1/SNR. Since this is
                        % done after despreading but before Viterbi
                        % decoder, we compensate for the code rate.
                        B = (P.NumberTxAntennas * P.RakeFingers) / (SNRlin * P.ConvRate) * eye(size(A));
                        G = (A + B) \ H_ii';
                        sTilde = G * VirtualAntennas;
                        
                    case 'SIC'
                        H_User = squeeze(H(:,:,ii));
                        yi = VirtualAntennas;
                        Hi = H_User;
                        % BPSK modulation
                        Constellations = [-1 1];   
                        for kk = 1:P.NumberTxAntennas
                            Gi = (Hi' * Hi) \ Hi';
                            g1i_star = Gi(1,:);
                            temp = g1i_star*yi;
                            % here we actually produce sHats and not sTilde
                            for bb = 1:NumOfEncBits/Users
                                [~, closestIndex] = min(abs(temp(bb) - Constellations));
                                sTilde(kk,bb) = Constellations(closestIndex);
                            end
                            yi = yi - Hi(:,1) * sTilde(kk);
                            Hi(:,1) = [];
                        end 
                    
                    otherwise
                        error('MIMO Detector not supported');
                end
                
                % since all antennas sent the same data, we can average the
                % estimates of each Tx antenna
                mrc = mean(sTilde, 1);
                
                % Decoding the bits: soft Viterbi decoder
                decodedBits  = step(decoder, real(mrc).');
                
                % Eliminating convolution tails
                RxBits(:,ii) = decodedBits(1:P.BitsPerUser).';
            end

            % Flatten the bit vectors for BER count
            RxBits      = reshape(RxBits, NumOfBits, 1);
            % Error count
            errors      = sum(RxBits ~= Bits);
            % Add to the errors found in the previous frames
            Results(ss) = Results(ss) + errors;

        end
    end
    
    % Calculate the BER
    BER = Results/(NumOfBits*P.NumberOfFrames);
end