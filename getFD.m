function [FD] = getFD(time_data, method)
% Input Parameters:
% Can different methods 'PS' = Power Spectra method
                    % 'Mono' = Monofractal method from Eke et al (2000)
                  % 'FDTool' = FD toolbox using Fractal Volatility
                  
TR = 1/250; % 1/freq (in seconds);
% method = 'PS';

% Finding length of time series:
TT = length(time_data);
if TT == 2^nextpow2(TT)
    powertwo = nextpow2(TT);
else
    powertwo = nextpow2(TT) - 1;
end
bins = powertwo - 2;
databins = 2^(bins+2);

% Fractal Dimension computation:
switch method
    case 'PS'
        % ***** Power spectrum method: *****
        sumcheck = sum(time_data);
        if sumcheck == 0,                               % No signal, set FD to noise level = 1.5            
            FD = 1.5;
            sumcheck = [];
           
        else
            time_data = detrend(time_data);
            Yps = fft(time_data, databins);             % FFT of temporal time points
            Pps = abs(Yps).^2;                          % Power = signal squared
            freq = (1/TR)*((1:(databins/2))/databins);
            LA = log10(Pps)'; LF = log10(freq);
            LA = LA(1,1:(databins/2));
            
            % Fit the power spectra:
            [fits, statsPS]=robustfit(LF,LA);
            
            beta = (-1)*fits(2,1);
            hurst = (beta+1)/2;                         % Assuming fGn not fBm
            p_val = statsPS.p(2,1);                     % P value of fit
            FD = 2 - hurst;
        end
        
    case 'Mono'
         % ***** Eke's Monofractal method: *****
        [hurst,Rsquare] = getHurst_2(time_data,TR);
        
        % Hurst(isnan(Hurst)) = 0;
        hurst(hurst==0) = 1/2;        
        FD = 2 - hurst;  
        
    case 'FDTool'
        % ***** Fractal volatility method: *****
        [FD,FDstd] = fractalvol(time_data);
        
        FD(isnan(FD)) = 0;
        FD(FD==0) = 1/2;
        
    otherwise
        display ('You did not enter a method - please enter one and try again!')
end