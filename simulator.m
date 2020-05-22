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
            case 'Bypass'
                H = ones(Users, P.NumberRxAntennas*P.ChannelLength, P.NumberTxAntennas);
                % himp = ones(Users,1);
            case 'AWGN'
                H = ones(Users, P.NumberRxAntennas*P.ChannelLength, P.NumberTxAntennas);
                % himp = ones(Users,1);
            case 'Multipath'
                H = sqrt(1/2) * (...
                      randn(Users, P.NumberRxAntennas*P.ChannelLength, P.NumberTxAntennas) +...
                      1i*randn(Users, P.NumberRxAntennas*P.ChannelLength, P.NumberTxAntennas)...
                    );
                % himp = sqrt(1/2)* ( randn(Users,P.ChannelLength) + 1i * randn(Users,P.ChannelLength) );
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
                    y = zeros(1,NumOfChipsPerUser+P.ChannelLength-1,Users);
                    for i = 1:Users
                        y(1,:,i) = conv(mwaveform(1,:,i),himp(i,:)) + noise(1,:,i); 
                    end
                otherwise
                    error('Channel not supported')
            end

            % Rake receiver and decoder
            RxBits = zeros(Users,NumOfBits/Users);
            
            for rr = 1:Users
                UserSequence = SpreadSequence(:,rr);
                
                fingers = zeros(P.ChannelLength,NumOfEncBits/Users);

                for i = 1:P.RakeFingers
                    data    =  y(1,i:i+NumOfChipsPerUser-1,rr)./PNSequence.'; 
                    rxvecs  = reshape(data,SeqLen,NumOfEncBits/Users);
                    fingers(i,:) = 1/SeqLen * UserSequence.' * rxvecs;
                end

                % Symbols for soft decoder
                mrc = (1/norm(himp(rr,:))) * conj(himp(rr,:)) * fingers;

                % Decoding the bits: soft Viterbi decoder
                DecodedBits = step(decoder, real(mrc).');
                
                % Eliminating convolution tails
                RxBits(rr,:) = DecodedBits(1:P.BitsPerUser);
                
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