function [Graw,f,Y,U] = myspa(y,u,Fs)
%% 
% Computes a raw (unsmoothed) estimate of the transfer function
% from FFTs of the (scalar) input/output data.
%
% Inputs
%    y = N-by-1 vector of output data
%    u = N-by-1 vector of input data
%    Fs = Sampling frequency in Hz  = 1/(sample time)
%
% Outputs
%    G = Nw-by-1 complex frequency response 
%    f = Nw-by-1 frequency vector (Hz)
%    Y = Nw-by-1 FFT of output.  Only response for f>0 is returned.
%    U = Nw-by-1 FFT of input.  Only response for f>0 is returned.


%% Compute FFTs of input and output
% This is a two-sided spectrum, i.e. positive & negative frequency
N = numel(u);
NFFT = 2^nextpow2(N); % Next power of 2 from length of y
f = Fs/2*linspace(0,1,NFFT/2+1)';       % Frequency vector, Hz
w = f*2*pi;                             % Frequency vector, rad/sec
Y2 = fft(y,NFFT)/N;                
U2 = fft(u,NFFT)/N;

% Positive frequency response
idx = 1:NFFT/2+1;
Y = Y2(idx);
U = U2(idx);

% % Single-sided amplitudes
% Ymag = 2*abs(Y(idx));
% Umag = 2*abs(U(idx));

%% Compute (raw) empirical transfer functions

% Raw frequency response from FFTs of input and output
Graw = Y(idx)./U(idx);
