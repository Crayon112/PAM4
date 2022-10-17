clc;
clear;
close all;

useRaiseCosineFilter = 0;
useFWHM = 0;

bitsPerSymbol = 2;
bitRate = bitsPerSymbol * 18e9;
DSPSampleRate = 90e9;
PRBS = 2^16;
nSamplesPerBit = 16;
nBits = bitsPerSymbol * PRBS;
nSamplesPerSymbol = 2 * nSamplesPerBit;
nSymbols = nBits / bitsPerSymbol;
symbolRate = bitRate / bitsPerSymbol;
numTimeSamples = nSamplesPerSymbol * nSymbols;
timeWindow = (1 / symbolRate) * nSymbols;
dt = timeWindow / numTimeSamples;
time = (0:numTimeSamples - 1) .* dt;
f = (1 / dt) .* (-numTimeSamples / 2:numTimeSamples / 2 - 1) ./ numTimeSamples;
w = 2 * pi / timeWindow * [0:numTimeSamples / 2, -numTimeSamples / 2 + 1: -1]';


load('symbol_pattern.mat');
symbol_pattern = [symbol_pattern, symbol_pattern];


S = zeros(1, length(symbol_pattern)/2);
for k = 1: length(symbol_pattern) / 2
    if symbol_pattern(2*k-1) == 0 && symbol_pattern(2*k) == 0
        S(k) = 0;
    elseif symbol_pattern(2*k-1) == 1 && symbol_pattern(2*k) == 0
        S(k) = 1;
    elseif symbol_pattern(2*k-1) == 1 && symbol_pattern(2*k) == 1
        S(k) = 2;
    elseif symbol_pattern(2*k-1) == 0 && symbol_pattern(2*k) == 1
        S(k) = 3;
    end
end

symbol_pattern = S;


if (useRaiseCosineFilter == 0)
    SAMPLES = zeros(1, numTimeSamples);
    div = numTimeSamples / length(symbol_pattern);
    for index = 1:length(symbol_pattern)
        SAMPLES((index - 1)*div+1:index*div) = ones(1, div) ...
            .* symbol_pattern(index);
    end
    Idrive = SAMPLES ./ max(SAMPLES);
else
    data = resample(symbol_pattern, 2, 1, 8);
    Idrive = interpft(data, numTimeSamples);
end

Idrive = Idrive - mean(Idrive);
eyediagram(Idrive(1: 65536), 32);
SMP = Idrive(1: 16: end);
eyediagram(SMP(1: 16384), 16);

N = length(SMP);
ts = 1 / DSPSampleRate;
fs = 1 / ts;
T = N * ts;

if (useFWHM > 0)
    F_0 = useFWHM / (2 * sqrt(log(2)));
    f_VEGA = [0.0001 : (N / 2 - 1) + 0.0001, (-N / 2) + 0.0001 : (-1) + 0.0001] * fs / N;
    H_f = sin(pi.*f_VEGA./fs) ./ (pi .* f_VEGA ./ fs);
    v_A_VEGA = real(ifft(fft(SMP)./H_f)); % Apply pre-emphasis
else
    v_A_VEGA = SMP;
end


figure()
plot(v_A_VEGA(1: end), '.b')

v_A_VEGA = v_A_VEGA(1: 2: end);
figure();
plot(v_A_VEGA, '.r');

writematrix(v_A_VEGA.', 'PAM4.csv');
