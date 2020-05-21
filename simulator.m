% EPFL - Advanced Wireless Receivers
% Final Project:
% CDMA IS95 standard, system simulator
% Brian Odermatt, Francesco Gallo
% May 2020

function BER = simulator(P)

    if P.CDMAUsers > P.HadLen
       error('More users than available sequences')
    end
    
    % Initial values and constants
    ConstraintLength = P.ConstrLen;
    Users = P.CDMAUsers;
    
    HadamardMatrix = hadamard(P.HadLen)/sqrt(P.HadLen);
    SpreadSequence = HadamardMatrix;
    SeqLen         = P.HadLen;
    
    NumOfBits            = (P.BitsPerUser * P.Modulation)* Users;                 % per Frame
    NumOfEncBits         = (NumOfBits + (ConstraintLength-1)*Users)/P.ConvRate ;  % per Frame    
    NumOfEncBitsPerUser  = (NumOfBits/Users + ConstraintLength-1)/P.ConvRate;     % per user, per Frame
    NumOfChipsPerUser    = NumOfEncBits * SeqLen / Users;                         % per user, per Frame
    
    % generating PN sequence
    PNSequence     = genPN(NumOfChipsPerUser);

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
            NumOfRXChipsPerUser = NumOfChipsPerUser+P.ChannelLength-1;
        otherwise
            NumOfRXChipsPerUser = NumOfChipsPerUser;
    end
    
    
    Results = zeros(1,length(P.SNRRange));


    for ii = 1:P.NumberOfFrames

        disp(['Simulating: ' sprintf('%3.f', ii/P.NumberOfFrames*100) '%'])

        Bits = randi([0 1], NumOfBits/Users, Users); % Random Data

        % Convolutional encoding: rate 1/2
        EncBits = zeros(NumOfEncBits/Users, Users);
        for i = 1:Users
            EncBits(:,i) = step(encoder, Bits(:,i));
        end
        EncBits = EncBits.';

        % Modulation
        switch P.Modulation % Modulate Symbols
            case 1 % BPSK
                symbols = -(2*EncBits - 1);
            otherwise
                error('Modulation not supported')
        end

        % Orthogonal spreading
        txsymbols = SpreadSequence(:,1:Users) * symbols;

        % apply PN sequence
        waveform = txsymbols(:).*PNSequence;

        % reshape to add multi RX antenna suppport
        waveform  = reshape(waveform,1,NumOfChipsPerUser);
        mwaveform = repmat(waveform,[1 1 Users]);

        % Channel
        switch P.ChannelType
            case 'Bypass'
                himp = ones(Users,1);
            case 'AWGN'
                himp = ones(Users,1);
            case 'Multipath'
                himp = sqrt(1/2)* ( randn(Users,P.ChannelLength) + 1i * randn(Users,P.ChannelLength) );
            otherwise
                error('Channel not supported')
        end

        %%%
        % Simulation
        snoise = randn(1,NumOfRXChipsPerUser,Users) + 1i * randn(1,NumOfRXChipsPerUser,Users);

        % SNR Range
        for ss = 1:length(P.SNRRange)
            SNRdb  = P.SNRRange(ss);
            SNRlin = 10^(SNRdb/10);
            % normalize noise according to SNR and 
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

            % Receiver
            switch P.ReceiverType
                case 'Rake'

                    RxBits = zeros(Users,NumOfBits/Users);

                    for rr = 1:Users
                        UserSequence = SpreadSequence(:,rr);
                        
                        FrameLength = NumOfChipsPerUser;

                        fingers = zeros(P.ChannelLength,NumOfEncBitsPerUser);

                        for i = 1:P.ChannelLength
                            data    =  y(1,i:i+FrameLength-1,rr)./PNSequence.'; 
                            rxvecs  = reshape(data,SeqLen,NumOfEncBitsPerUser);
                            fingers(i,:) = 1/SeqLen * UserSequence.' * rxvecs;
                        end
                        
                        % Symbols for soft decoder
                        mrc = (1/norm(himp(rr,:))) * conj(himp(rr,:)) * fingers;
                        
                        % Decoding the bits
                        decodedBitsForUser = step(decoder, real(mrc).');
                        RxBits(rr,:) = decodedBitsForUser(1:P.BitsPerUser);
                    end

                otherwise
                    error('Receiver not supported')
            end
            
            % Flatten the bit vectors for BER count
            Bits    = reshape(Bits, NumOfBits, 1);
            RxBits  = reshape(RxBits.', NumOfBits, 1);

            % BER count
            errors =  sum(RxBits ~= Bits);
            Results(ss) = Results(ss) + errors;

        end
    end

    BER = Results/(NumOfBits*P.NumberOfFrames);
end