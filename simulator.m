% EPFL Advanced Wireless Receivers
% Project IS95, Spring 2020
% Francesco Gallo, Brian Odermatt


function BER = simulator(P)

    %% Setup parameters
    if P.CDMAUsers > P.HamLen
       disp('WARNING: More user then sequences');
       BER = -1;
       return;
    end
    
    
    RX = P.CDMAUsers;
    SeqLen = P.HamLen;
    NumberOfInformationBits = P.NumberOfSymbols;                    % per user
    NumberOfEncodedBits = P.Modulation * NumberOfInformationBits;   % per user
    NumberOfModulatedBits = NumberOfEncodedBits * P.HamLen;
    switch P.ChannelType
        case 'Multipath', NumberOfChipsRX = NumberOfModulatedBits + P.ChannelLength - 1;
        otherwise, NumberOfChipsRX = NumberOfModulatedBits;
    end
    
    % Convolutional encoder
    encoder = comm.ConvolutionalEncoder(...
        'TerminationMethod', 'Continuous',...
        'TrellisStructure', poly2trellis(9, [753 561])...
    );
    
    % Orthogonal modulation
    WalshMatrix = (-hadamard(64) + 1)/2;
    WalshFunctions = WalshMatrix(1:RX,:);
    
    % PN spreading polynomials
    % in-phase (I)
    iPwrs = [15, 13, 9, 8, 7, 5, 0]';
    iGen = zeros(16, 1);
    iGen(16-iPwrs) = ones(size(iPwrs));     % Gi = [ 1 0 1 0 0 0 1 1 1 0 1 0 0 0 0 1]';
    iState = [zeros(length(iGen)-1, 1); 1]; % Initial State
    
    % quadrature-phase (Q)
    qPwrs = [15, 12, 11, 10, 6, 5, 4, 3, 0]';
    qGen = zeros(16, 1);
    qGen(16-qPwrs) = ones(size(qPwrs));     % Gq = [ 1 0 0 1 1 1 0 0 0 1 1 1 1 0 0 1]';
    qState = [zeros(length(qGen)-1, 1); 1]; % Initial State

    
    % Convolutional decoder
    decoder = comm.ViterbiDecoder(...
        'TerminationMethod', 'Continuous',...
        'TracebackDepth', 5*9,...
        'TrellisStructure', poly2trellis(9, [753 561])...
    );

    
    

    for ii = 1:P.NumberOfFrames
        ii

        %% Transmitter
        % Generate random information bits for the user
        informationBits = randi([0, 1], RX, P.NumberOfSymbols);

        % Convolutional Encoding
        encodedBits = zeros(NumberOfEncodedBits, RX);
        for i = 1:RX
            encodedBits(:,i) = encoder(informationBits(i,:).');
        end
        encodedBits = encodedBits.';

        % Orthogonal Modulation
        modulatedBits = reshape(WalshFunctions.'*encodedBits, NumberOfModulatedBits, 1);

        % PN spreading (quadrature spread length 2^15)
        [iPN, iState] = lsfrPN(iGen, iState, NumberOfModulatedBits);
        [qPN, qState] = lsfrPN(qGen, qState, NumberOfModulatedBits);

        % modulation
        switch P.Modulation
            case 2  % QPSK
                pnModulated = sign(iPN-1/2) + 1i*sign(qPN-1/2);
            otherwise, disp('Modulation not supported')
        end
        
        wave = modulatedBits(:) .* pnModulated;

        % reshape to add multi RX antenna support (from 1 col-vector to
        % multiple copies of a row vector)
        wave = reshape(wave, 1, NumberOfModulatedBits);
        mwave = repmat(wave, [1 1 RX]);


        %% Channel
        switch P.ChannelType
            case 'AWGN', himp = ones(RX,1);
            case 'Multipath', himp = sqrt(1/2)*(randn(RX,P.ChannelLength) + 1i*randn(RX,P.ChannelLength));
            otherwise, disp('Channel not supported')
        end
        % random complex normal noise ~N(0,1)
        snoise = randn(1,NumberOfChipsRX,RX) + 1i*randn(1,NumberOfChipsRX,RX);

        % SNR Range
        for ss = 1:length(P.SNRRange)
            SNRdb  = P.SNRRange(ss);
            SNRlin = 10^(SNRdb/10);
            noise  = 1/sqrt(2*SeqLen*SNRlin)*snoise;

            % Channel: add noise and multipath
            switch P.ChannelType
                case 'Bypass', y = mwave;
                case 'AWGN', y = mwave + noise;
                case 'Multipath'     
                    y = zeros(1, NumberOfChipsRX, RX);
                    for i = 1:RX
                        y(1,:,i) = conv(mwave(1,:,i),himp(i,:)) + noise(1,:,i); 
                    end
                otherwise, disp('Channel not supported')
            end     
         end
    end

end