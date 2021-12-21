# Pressure Signal from RA6 Reactor
The crow data were analyzed in the frequency domain using the power spectral density (PSD). For this, the Welch method has been used, given the sampling time Ts= 0,01s and taking a Nfft = 1204 to perform the calculation of the discrete Fourier transform. In this way, a frequency resolution of: F = 1/(Nfft*Ts) = 0,09Hz.
The Welch method allows to estimate the PSD of a large amount of data to be calculated. The procedure consists of calculating a certain amount of PSDs and then averaging them, where each PSD is not independent of each other, since they are calculated with overlapping data.
At the same time, distribution and density probability have been calculated:
scipy and numpy libraries were used.
