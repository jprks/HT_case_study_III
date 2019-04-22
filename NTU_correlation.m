function x_ntu = NTU_correlation(x)
x_ntu = 1 - exp(x^0.22*(exp(-x^0.78)-1));
end