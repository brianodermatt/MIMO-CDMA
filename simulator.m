% Wireless Receivers II - Assignment 2:
%
% Direct Sequence Spread Spectrum Simulation Framework
%
% Telecommunications Circuits Laboratory
% EPFL


function BER = simulator(P)

    if P.CDMAUsers > P.HamLen
       disp('WARNING: More user then sequences');
       BER = -1;
       return;
    end
    RX = P.CDMAUsers;
    
    % Generate the spreading sequences
    HadamardMatrix = hadamard(P.HamLen)/sqrt(P.HamLen);            
    SpreadSequence = HadamardMatrix;
    
    SeqLen         = P.HamLen;
    
    NumberOfChips  = P.NumberOfSymbols*P.Modulation*SeqLen; % per Frame

    PNSequence     = genbarker(NumberOfChips); % -(2*step(GS)-1);
    
    NumberOfBits   = P.NumberOfSymbols*P.Modulation*RX; % per Frame
    
    encoder = comm.ConvolutionalEncoder(...
        'TerminationMethod', 'Continuous',...
        'TrellisStructure', poly2trellis(9, [753 561])...
    );
    decoder = comm.ViterbiDecoder(...
        'TerminationMethod', 'Continuous',...
        'TracebackDepth', 5*9,...
        'TrellisStructure', poly2trellis(9, [753 561])...
    );
    
    % Channel
    switch P.ChannelType
        case 'Multipath',
            NumberOfChipsRX   = NumberOfChips+P.ChannelLength-1;
        otherwise,
            NumberOfChipsRX = NumberOfChips;
    end

Results = zeros(1,length(P.SNRRange));

for ii = 1:P.NumberOfFrames
    ii
    
    %% Encoder
    % Generate random information bits for the user
    informationBits = randi([0, 1], RX, NumberOfBits/RX);
    
    % Convolutional Encoding
    for i = 1:RX
        encodedBits(:,i) = encoder(informationBits(i,:).');
    end
    encodedBits = encodedBits.';
    
    % Orthogonal Modulation
    WalshMatrix = (-hadamard(64) + 1)/2;
    WalshFunctions = WalshMatrix(1:RX,:);
    modulatedBits = WalshFunctions.'*encodedBits;
    
    % PN spreading (quadrature spread length 2^15)
    
    %% Old code
    bits = randi([0 1],1,NumberOfBits); % Random Data
 
    % Modulation
    switch P.Modulation % Modulate Symbols
        case 1, % BPSK
            symbols = -(2*bits - 1);
        otherwise,
            disp('Modulation not supported')
    end
    

    % distribute symbols on users
    SymUsers = reshape(symbols,RX,NumberOfBits/RX);
    
    % Convolution Enconding
    
    % multiply hadamard
    txsymbols = SpreadSequence(:,1:RX) * SymUsers;
        
    % apply Barker code
    waveform = txsymbols(:).*PNSequence;

    % reshape to add multi RX antenna suppport
    waveform  = reshape(waveform,1,NumberOfChips);
    mwaveform = repmat(waveform,[1 1 RX]);
    
    % Channel
    switch P.ChannelType
        case 'AWGN',
            himp = ones(RX,1);
        case 'Multipath',
            himp = sqrt(1/2)* ( randn(RX,P.ChannelLength) + 1i * randn(RX,P.ChannelLength) );
%             himp = (ones(RX,1) * sqrt(P.PDP)) .* himp;
        otherwise,
            disp('Channel not supported')
    end
    
    %%%
    % Simulation
    snoise = ( randn(1,NumberOfChipsRX,RX) + 1i* randn(1,NumberOfChipsRX,RX) );
    
    % SNR Range
    for ss = 1:length(P.SNRRange)
        SNRdb  = P.SNRRange(ss);
        SNRlin = 10^(SNRdb/10);
        noise  = 1/sqrt(2*SeqLen*SNRlin) *snoise;
        
        % Channel
        switch P.ChannelType
            case 'AWGN',
                y = mwaveform + noise;
            case 'Multipath'     
                y = zeros(1,NumberOfChips+P.ChannelLength-1,RX);
                for i = 1:RX
                    y(1,:,i) = conv(mwaveform(1,:,i),himp(i,:)) + noise(1,:,i); 
                end
            otherwise,
                disp('Channel not supported')
        end
        
 
        % Receiver
        switch P.ReceiverType
            case 'Rake',
                rxbitsuser = zeros(RX,NumberOfBits/RX);
                
                for rr=1:RX
                    UserSequence = SpreadSequence(:,rr);
                    
                    
                    FrameLength = P.NumberOfSymbols * SeqLen;
                    
                    fingers = zeros(P.ChannelLength,P.NumberOfSymbols);
                
                    for i=1:P.ChannelLength
                        data    =  y(1,i:i+FrameLength-1,rr).*PNSequence.'; 
                        rxvecs  = reshape(data,SeqLen,P.NumberOfSymbols);
                        fingers(i,:) = 1/SeqLen * UserSequence.' * rxvecs;
                    end
                    mrc = (1/norm(himp(rr,:))) * conj(himp(rr,:)) * fingers;
                    rxbitsuser(rr,:) = real(mrc) < 0;
                end
                rxbits = reshape(rxbitsuser,1,NumberOfBits);
            otherwise,
                disp('Source Encoding not supported')
        end
        
        % BER count
        errors =  sum(rxbits ~= bits);
        
        Results(ss) = Results(ss) + errors;
        
    end
end

BER = Results/(NumberOfBits*P.NumberOfFrames);
end

function seq = genbarker(len)
    BarkerSeq = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];

    factor = ceil(len/length(BarkerSeq));
    b = repmat(BarkerSeq,1,factor);
    b = BarkerSeq.'*ones(1,factor);
    seq = b(1:len).';
end