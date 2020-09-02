[samples, h, p] = BoltzmannSampling(10000, 10000);

function [samples, h, p] =  BoltzmannSampling(M,N)

%%% Create empty vectors for the output %%%
samples = zeros(M,N);
h = zeros(N,1);
p = zeros(N,1);

%%% Set value for A %%%
A = (4*exp(-(1/2)))/sqrt(pi);


%%% Set estimate for number of samples needed to end with M 
%%% accepted values. We drop samples if we overshoot M and we add samples if
%%% we undershoot M. This will speed up the program by decreasing the time
%%% spent in the 'while loop' below. %%% 
est = round(M*A);

%%% Begin 'for loop' for the N trials %%%
for i = 1:N
    
    %%% Sample random variables from Wiebull distribution using the inverse
    %%% CDF of the Wiebull distribution. %%%
    X = 2.*sqrt(-log(1-rand(est,1)));
    Z = rand(est,1);
    
    %%% Creates a vector of logical values of 1 and 0, where 1 is the 
    %%% acceptance of X_i as Maxwellian distributed
    Y = (Z.*A.*(X./2).*exp(-(X./2).^2) <= sqrt(2/pi).*(X.^2).*exp(-(X.^2)./2));
    
    %%% 'While loop to add samples until we have M accepted values
    %%% for Maxwellian distributed X_i. %%%
    while sum(Y) < M
        x = 2.*sqrt(-log(1-rand(1)));
        z = rand(1);
        y = (z*A*(x/2)*exp(-(x/2)^2) <= sqrt(2/pi)*(x^2)*exp(-(x^2)/2));
        Y = [Y ; y ];
        X = [X ; x];
    end
    
    %%% Acc are the accepted values %%%
    Acc = nonzeros(X.*Y);
    
    %%% Remove extra accepted values so there are M = 10,000 samples. %%%
    Acc(M+1:end) = [];
    
    %%% Create Maxwellian samples generated from standard normals %%%
    Q = sqrt(normrnd(0,1, [M,1]).^2 + normrnd(0,1, [M,1]).^2 + normrnd(0,1,[M,1]).^2 );
    
    %%% Output values from KS-test of the two Maxwellian distributed samples. %%%
    [h(i), p(i)] = kstest2(Acc, Q); 
    
    %%% Store output values of the Maxwellian distributed samples generated
    %%% using rejection sampling. %%%
    samples(:,i) = Acc;
    
end
end