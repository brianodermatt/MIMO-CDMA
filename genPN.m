function seq = genPN(len)
    
    % Generating PN sequence following CDMA IS-95 standard

    % PN spreading polynomials
    % in-phase (I)
    iPwrs = [15; 13; 9; 8; 7; 5; 0];
    iGen = zeros(16, 1);
    iGen(16-iPwrs) = ones(size(iPwrs));     % Gi = [ 1 0 1 0 0 0 1 1 1 0 1 0 0 0 0 1]';
    iState = [zeros(length(iGen)-1, 1); 1]; % Initial State
    % quadrature-phase (Q)
    qPwrs = [15; 12; 11; 10; 6; 5; 4; 3; 0];
    qGen = zeros(16, 1);
    qGen(16-qPwrs) = ones(size(qPwrs));     % Gq = [ 1 0 0 1 1 1 0 0 0 1 1 1 1 0 0 1]';
    qState = [zeros(length(qGen)-1, 1); 1]; % Initial State
    % PN spreading (quadrature spread length 2^15)
    [iPN, ~] = lsfrPN(iGen, iState, len);
    [qPN, ~] = lsfrPN(qGen, qState, len);
    
    % IQ modulated PN sequence
    seq = - sign(iPN-1/2) - 1i*sign(qPN-1/2);
    
end

