function [pn, state] = lsfrPN(generator,initialState,N)
%lsfrPN
%   [spread, newState] = lsfrPN(generator,initialState,length)
%   Generates Pseudo Noise sequence using a linear feedback register
%   
%   Inputs:
%   generator - generation polynom: [1 0 1 0 1 ...]'
%   initialState - initial state of the Linear Shift Feedback Register
%   N - the length of the sequence to generate
%
%   Outputs:
%   pn - the generated PN sequence
%   state - the new state of the LSFR

    L = length(generator);
    state = initialState;
    
    pn = zeros(N,1);
    for i = 1:N
       pn(i) = state(L);
       state = xor(generator * state(L), state);
       state = [state(L); state(1:L-1)];
    end

end

