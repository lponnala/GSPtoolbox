% DETECTF: "Detect f=1/3"
% Tests for the presence of a sinusoid with frequency 1/3
% Usage: H = detectf(signal,alpha,display) 
% alpha = significance level (0.05 typically)
% (display=1): plot periodogram if sinusoid is detected
% If H=1, there exists a sinusoid of frequency 1/3
% If H=0, it does not exist

function H = detectf(signal,alpha,display)

% Compute periodogram
Y = fft(signal);
N = length(Y);
% If signal length is not a multiple of 3, then f=1/3 will not be a 
% Fourier frequency!
if rem(N,3)~=0
    error('\nSignal length is not a multiple of 3!');
end
I = (1/N)*(abs(Y).^2);

% Assign discrete frequencies
if rem(N,2)==0
    k = [0, (1:1:floor((N-1)/2)), (N/2), (floor((N-1)/2):-1:1)];
elseif rem(N,2)==1
    k = [0, (1:1:(N-1)/2), ((N-1)/2:-1:1)];
end

I_0 = I(find(k==0));
I_k = sum(I(find(k==N/3)))/2;

F = ((N-3)*I_k)/((signal*signal')-I_0-(2*I_k));
Fc = finv(1-alpha,2,N-3);
if F>Fc
    H = 1;
else
    H = 0;
end

if (H==1)&(display==1) % Change this to (H==0) if you want the periodogram to be displayed when f=1/3 is not detected
    if rem(N,2)==0
        ind = 2:1+(N/2); pgval = I(ind); freq = (1:(N/2))/N;
    elseif rem(N,2)==1
        ind = 2:1+((N-1)/2); pgval = I(ind); freq = (1:((N-1)/2))/N;
    end    
    plot(freq,pgval), grid on
    xlabel('cycles/base')
    title('Periodogram')
    pause
end

% % ----- OLD CODE -----
% if (H==1)&(display==1)
%     % Plot periodogram
%     Y(1) = [];
%     power = (1/N)*abs(Y(1:floor(N/2))).^2;
%     nyquist = 1/2;
%     freq = (1:N/2)/(N/2)*nyquist;
%     plot(freq,power), grid on
%     xlabel('cycles/base')
%     title('Periodogram')
%     pause(1)
% end