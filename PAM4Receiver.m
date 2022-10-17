clear;
close all;

global h_xx;
global R;
global Strat_P;


Strat_P = 1;
Syn_N = 1;

R_X = readmatrix('RX_X.csv');
inputSignal = R_X / mean(abs(R_X));

Resolution = 6;
inputSignal = quant(inputSignal, max(max(inputSignal))/2^Resolution);

Nor_Fac = (3 + 2 + 1) / 4;
inputSignalLevel = inputSignal / mean(abs(inputSignal)) * Nor_Fac;
X = fftshift(abs(fft(inputSignalLevel)).^2);

load('Traning_symbols.mat');
TS1 = Traning_symbols;
TS = TS1(1:end);

RX_Syn = inputSignalLevel;
RX_Syn = RX_Syn / mean(abs(RX_Syn)) * Nor_Fac;
RX_Syn = RX_Syn * 2 - 3;

TS = TS * 2 - 3;
R = [0, 1, 2, 3];
X_TS1 = TS;

for tap_number = 31
    choice = 1;
    h_xx(tap_number) = 0;
    h_xx(tap_number-fix(tap_number/2)) = 1;

    D_EQ = 3;
    step = 0.0005;
    for i = 1:30
        [out_X] = SP_EQ_DD(step, tap_number, RX_Syn, D_EQ, X_TS1);
    end
    Power_x = out_X(2:2:end);
    IQX = Power_x;
    [BER] = decision_BRR_counting(Power_x, X_TS1);
end


function [out_X] = SP_EQ_DD(step, tap_number, RX, D_EQ, X_TS1)
global h_xx;
global R;
global Strat_P;

X(tap_number) = 0;

piece = 500;
save_initial = 1;
id_E = 2;
choice = 1;

vary_error_x = zeros(1, length(RX)-piece);
vary_error_x(1) = 1;

X(save_initial) = RX(1);


if D_EQ == 3 %%  training sequence  (TMA)
    DX_Power = X_TS1;
    delay = ceil((tap_number - 1)/4);
    DX_Power = [zeros(1, delay+Strat_P), DX_Power];
end

X_all = RX;

for yyy = 1:30
    X(save_initial) = RX(1);
    out_X = zeros(1, length(X_all)-piece);
    error_x_T = zeros(1, length(X_all)-piece);
    error_x = zeros(1, length(X_all)-piece);
    index_X = zeros(1, length(X_all)-piece);

    for i = 1:piece:(length(X_all) - piece)
        for jj = 0:1:(piece - 1)
            out_X(i+jj) = X * h_xx.';
            if mod((i + jj), 2) == 0
                index_sym = (i + jj) / 2;
                if D_EQ == 1 % CMA
                    error_x(i+jj) = 1 - (abs(out_X(i+jj)));
                elseif D_EQ == 2 % MCMA
                    error_x_T(1) = R(1) - ((out_X(i+jj)));
                    error_x_T(2) = R(2) - ((out_X(i+jj)));
                    error_x_T(3) = R(3) - ((out_X(i+jj)));
                    error_x_T(4) = R(4) - ((out_X(i+jj)));
                    index_X(i+jj) = find(abs(error_x_T) == min(abs(error_x_T)));
                    error_x(i+jj) = error_x_T(index_X(i+jj));
                elseif D_EQ == 3 % TMA
                    error_x(i+jj) = DX_Power(index_sym) - (out_X(i+jj));
                end
            else
                error_x(i+jj) = 0;
            end

            % 抽头数更新
            h_xx = h_xx + step * error_x(i+jj) .* X;
            for save_index = 1:1:(tap_number - 1)
                X(tap_number-save_index+1) = X(tap_number-save_index);
            end
            X(save_initial) = X_all(i+jj+1);
        end

        vary_error_x(id_E) = mean(abs(error_x(i+1:1:i+jj)).^2);
        if choice == 1
            if vary_error_x(id_E) <= 0.05
                choice = 2;
                step = 0.00051;
            end
        end
        id_E = id_E + 1;
    end
end
end


function [BER] = decision_BRR_counting(Power_x, X_TS1)
TH_X2 = mean(Power_x);
NX1 = Power_x > TH_X2;
TH_X3 = mean(Power_x(NX1));
NX2 = Power_x < TH_X2;
TH_X1 = mean(Power_x(NX2));

Des_S(Power_x <= TH_X1) = -3;
Des_S(Power_x > TH_X1 & Power_x <= TH_X2) = -1;
Des_S(Power_x > TH_X2 & Power_x <= TH_X3) = 1;
Des_S(Power_x > TH_X3) = 3;

Syn_TS_X = X_TS1(1:2000);
N_s = 2000;
XP = zeros(1, 1000);
XS = zeros(1, 1000);
XM = zeros(1, 1000);
for d = 1:1000
    XP(d) = sum(conj(Des_S(d:d+N_s-1)).*Syn_TS_X);
    XS(d) = sqrt(sum(abs(Des_S(d:d+N_s-1)).^2)*sum(abs(Syn_TS_X).^2));
    XM(d) = abs(XP(d)/XS(d))^2;
end

Start_p_X = find(XM == max(XM));


I_X1 = [Des_S(Start_p_X:end), Des_S(1:Start_p_X-1)];

[error, error_bit_low, error_bit_middle, error_bit_high, error_Sam] = BRR_counting(I_X1, X_TS1);
figure()
subplot(3, 1, 1)
plot(error_bit_low, 'o-g');
subplot(3, 1, 2)
plot(error_bit_middle, 'o-b');
subplot(3, 1, 3)
plot(error_bit_high, 'o-r');
Biterror = error;
Biterror1 = error_bit_low;
Biterror2 = error_bit_middle;
Biterror3 = error_bit_high;
BER = (sum(Biterror)) ./ (length(Biterror) * 2);
BER_low = (sum(Biterror1)) ./ (length(Biterror1));
BER_middle = sum(Biterror2) ./ (length(Biterror2));
BER_high = sum(Biterror3) ./ (length(Biterror3));

error_Sam_num = sum(error_Sam);
fprintf(['\n  BER_low=', num2str(BER_low), ' , BER_middle=', num2str(BER_middle), ' , BER_high=', num2str(BER_high)]);
fprintf(['\n error_Sam=', num2str(error_Sam_num)]);
fprintf(['\n BER=', num2str(BER)]);
end

function [error, error_bit_low, error_bit_middle, error_bit_high, error_Sam] = BRR_counting(I_X, X_TS)

error_Sam = zeros(1, length(I_X));
for k = 1:length(I_X)
    if I_X(k) ~= X_TS(k)
        error_Sam(k) = 1;
    else
        error_Sam(k) = 0;
    end
end

Decision_x = zeros(length(I_X), 2);
D_X = zeros(length(I_X), 2);
for i = 1:length(I_X)
    if I_X(i) == 3
        Decision_x(i, :) = [0, 1];
    elseif I_X(i) == 1
        Decision_x(i, :) = [1, 1];
    elseif I_X(i) == -3
        Decision_x(i, :) = [0, 0];
    elseif I_X(i) == -1
        Decision_x(i, :) = [1, 0];
    end

    if X_TS(i) == 3
        D_X(i, :) = [0, 1];
    elseif X_TS(i) == 1
        D_X(i, :) = [1, 1];
    elseif X_TS(i) == -3
        D_X(i, :) = [0, 0];
    elseif X_TS(i) == -1
        D_X(i, :) = [1, 0];
    end
end

error = zeros(1, length(I_X));
for j = 1:length(I_X)
    if Decision_x(j) ~= D_X(j)
        error(j) = 1;
    else
        error(j) = 0;
    end
end

error_bit_low = zeros(1, length(I_X));
error_bit_middle = zeros(1, length(I_X));
error_bit_high = zeros(1, length(I_X));

for i = 1:length(Decision_x(:, 1))
    if Decision_x(i, 1) ~= D_X(i, 1)
        if Decision_x(i, 2) == 0
            error_bit_low(i) = 1;
        elseif Decision_x(i, 2) == 1
            error_bit_high(i) = 1;
        end
    end

    if Decision_x(i, 2) ~= D_X(i, 2)
        error_bit_middle(i) = 1;
    else
        error_bit_middle(i) = 0;
    end
end

end
