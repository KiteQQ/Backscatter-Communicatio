# Backscatter-Communicatio

clc;
clear;

N=5000;  %N*L(The number of RF source transmitted symbols)
alpha = 0.001;   %the square of reflection parameter

M = 2; %The number of receive antennas
j = sqrt(-1);
ALP = 1/sqrt(2)* [1+1j,1-1j ,-1+1j ,-1-1j]; % All possible transmitted symbols from RF source
ALC = [0,1];  %All possible transmitted symbols from BD
Q = 4;  %The RF source modulation order
fs = 200;   %spreading gain

SNR_dB = [16:2:34]; % SNR in dB of RF source
channel = 10000000;   %The number of Monte Carlo realizations
varNoise = 1;  %noise power

for kChannel = 1:channel
    h1 = 1/sqrt(2) * (normrnd(0,1,M,1) + j*normrnd(0,1,M,1));       %channel gain = 1 channel from RF source to reader
    h2 = 1/sqrt(2) * (normrnd(0,1,M,1) + j*normrnd(0,1,M,1));  %channel gain = 1 channel from BD source to reader
    g = 1/sqrt(2) * (normrnd(0,1,1,1) + j*normrnd(0,1,1,1));   % channel from RF source to BD

    for kSNR = 1:size(SNR_dB,2)
        
        %%  SIGNAL MODEL
        P(kSNR) = 10^(SNR_dB(kSNR)/10)*varNoise;
        s = 1/sqrt(2) * ([sign(rand(N,1)-0.5)] + j* [sign(rand(N,1)-0.5)]); %RF source transmitted signal
        c =[-1;1;sign(rand(N/(fs)-2,1)-0.5)];  %BD transmitted signal
        c(c==-1)=0;
        c_sample = reshape((repmat(c,1,fs)).',N,1);   %
        noise = 1/sqrt(2) * (normrnd(0,sqrt(varNoise),M,(N)) + j*normrnd(0,sqrt(varNoise),M,(N)));  %noise
        y =  sqrt(P(kSNR))*h1 * s.'  +  sqrt(P(kSNR)* alpha * g  )* h2 * ((s.') .*( c_sample.'))+ [noise];
        
        %% channel estimation
        z1 = abs(y(:,1:fs)).^2.*exp(j*Q*angle(y(:,1:fs)));  % the phase statistics  corresponding to c=0
        theta1 = (angle( -sum(z1,2)/fs))/(Q); % the phase corresponding to c=0 using M-power
        R1 = sqrt((sum(abs(y(:,1:fs)).^2,2))/(fs)-1); % The magnitude estimation corresponding to c=0

        C1=0;
        for k = 1:fs
            C1 = C1+ y(:,k)*y(:,k)';  % The sample covariance corresponding to c=0
        end
          [e1,v1] = eig(C1/fs-diag(ones(M,1))); % The eigenvalue and eigenvector corresponding to c=0
         [a,b] = max(diag(v1)); % The maximum eigenvalue corresponding to c=0
         h1_mod =sqrt(a)* e1(:,b); % estimated channel without phase correction corresponding to c=0
         
           [a,b] =  max(R1);  % The maximum power for each antenna corresponding to c=0
           
         h1_mod_modified = h1_mod(b) / exp(j*theta1(b))/abs(h1_mod(b) ); % Phase correction corresponding to c=0
         h1_hat = h1_mod/h1_mod_modified;  % estimated channel with phase correction corresponding to c=0
       
        z2 = abs(y(:,fs+1:2*fs)).^Q.*exp(j*Q*angle(y(:,fs+1:2*fs)));  % the phase statistics  corresponding to c=1
        theta2 = (angle(-sum(z2,2)/fs))/(Q);  % the phase corresponding to c=1 usign M-power
        R2 = sqrt((sum(abs(y(:,fs+1:2*fs)).^2,2))/(fs)-1); % The magnitude estimation corresponding to c=1

        C2=0;
        for k = fs+1:2*fs
            C2 = C2+ y(:,k)*y(:,k)'; % The sample covariance corresponding to c=1
        end
        [e2,v2] = eig(C2/fs-diag(ones(M,1))); % The eigenvalue and eigenvector corresponding to c=1
        [a,b] = max(diag(v2));   % The maximum eigenvalue corresponding to c=1
        h2_mod = sqrt(a)*e2(:,b); % estimated channel without phase correction corresponding to c=1
        
        
         [a,b] =  max(R2); % The maximum power for each antenna corresponding to c=1
         h2_mod_modified = h2_mod(b) / exp(j*theta2(b))/abs(h2_mod(b) );    % Phase correction corresponding to c=1
         h2_hat =   h2_mod/h2_mod_modified;   % estimated channel with phase correction corresponding to c=1
       
%% P-SD detection
        reconstruction_label3 = [h1_hat*ALP, h2_hat*ALP];
        [px] = optimalDetectorforN50(y.',8,reconstruction_label3.',1); %LLRT detector
        px2 = sum(px(:,1:4),2);
        px3 = sum(px(:,5:8),2);
        c_est_power(1:2,:)=[0;1];
        for   k = 3:N/(fs)
            if sum(sum(log(px2(fs*(k-1)+1:fs*k,:)./px3(fs*(k-1)+1:fs*k,:))))>0
                c_est_power(k,:) = 0;
            else
                c_est_power(k,:) = 1;
            end
        end

        err_m_power_psd(kSNR,kChannel) =  1/(N/(fs)-2) * ( sum( c_est_power(3:end) ~= c(3:end))); %BER for P-SD

%% JCD
        max_iter_num =300; % The maximum iteration number
        threshold = 0.01;  % The threshold of convergence
        [c_est_power,ite] = M_Power_iteration(c_est_power,max_iter_num,N,fs,y,Q,M,ALM,ALP,threshold);  %JCD algorithm 
        ite_num(kSNR,kChannel) = ite; %The returned iteration number
        
        err_m_power(kSNR,kChannel) =  1/(N/(fs)-2) * ( sum( c_est_power(3:end) ~= c(3:end))); %BER for JCD
        
       
        
        %%   Perffect CSI
        reconstruction_label = [(sqrt(P(kSNR))*h1)*ALP, (sqrt(P(kSNR))*h1+sqrt(P(kSNR) * alpha * g )* h2)*ALP];
        [px6] = optimalDetectorforN50(y.',8,reconstruction_label.',1);
        px7 = sum(px6(:,1:4),2);
        px8 = sum(px6(:,5:8),2);
        for   k = 1:N/(fs)
            if sum(log(px7(fs*(k-1)+1:fs*k)./px8(fs*(k-1)+1:fs*k)))>0
                c_est_ml(k,:) = 0;
            else
                c_est_ml(k,:) = 1;
            end
        end
        err_ml(kSNR,kChannel) =  1/(N/(fs)-2) * ( sum( c_est_ml(3:end) ~= c(3:end)));  %BER for Optimal Detector
        
    end
end

%%
average_err_m_power_psd = sum(err_m_power_psd,2)/channel;
average_err_m_power = sum(err_m_power,2)/channel;
average_err_ml = sum(err_ml,2)/channel;
average_ite_num = sum(ite_num,2)/channel;

figure
semilogy(SNR_dB,average_err_m_power_psd ,'b*-', 'LineWidth',1)
hold on
semilogy(SNR_dB,average_err_m_power,'ro-', 'LineWidth',1)
hold on
semilogy(SNR_dB,average_err_ml,'k<-', 'LineWidth',1)
xlabel('SNR (dB)'), ylabel('BER')
grid on
legend('M-Power','Blind','ML')
title('QPSK N=50 K=4 -30')

figure
plot(SNR_dB,average_ite_num)





function   [c, ite]  =  M_Power_iteration(c,iter_number,N,fs,y,Q,M,ALP,threshold)
%% Algorithm for JCD

  h1_hat1 = 0; % Initianlized value
  h2_hat1 = 0;
for ite = 1:iter_number
    c_est_power = c;
    z3 = [];  % The initialized set corresponding to c=0
    z4 = [];   % The initialized  set corresponding to c=1
    for k = 1:N/fs
        if c_est_power(k,:) == 0
            z3 =[z3, y(:,(k-1)*fs+1:k*fs)];
        else
            z4 =[z4, y(:,(k-1)*fs+1:k*fs)];
        end
    end
    
    % After updating the data set corresponding tp c=0 and c=1 respectively, the process of JCD is the same with P-SD
    z5 = abs(z3).^2.*exp(j*Q*angle(z3));  
    theta1 = (angle( -sum(z5,2)/fs))/(Q);
    R1 = sqrt((sum(abs(z3).^2,2))/(size(z3,2))-1);
    
    C1=0;
    for k = 1:size(z3,2)
        C1 = C1+ z3(:,k)*z3(:,k)';
    end
    [e1,v1] = eig(C1/size(z3,2)-diag(ones(M,1)));
    [a,b] = max(diag(v1));
    h1_mod =sqrt(a)* e1(:,b);

    [a,b] =  max(R1);
    
    h1_mod_modified = h1_mod(b) / exp(j*theta1(b))/abs(h1_mod(b) );
    h1_hat = h1_mod/h1_mod_modified;
    
     z6 = abs(z4).^2.*exp(j*Q*angle(z4));
    theta2 = (angle(-sum(z6,2)/fs))/(Q);
    
    R2 = sqrt((sum(abs(z4).^2,2))/(fs)-0.8);
    C2=0;
    for k = 1:size(z4,2)
        C2 = C2+ z4(:,k)*z4(:,k)';
    end
    [e2,v2] = eig(C2/size(z4,2)-diag(ones(M,1)));
    [a,b] = max(diag(v2));
    h2_mod = sqrt(a)*e2(:,b);
    
    [a,b] =  max(R2);
    h2_mod_modified = h2_mod(b) / exp(j*theta2(b))/abs(h2_mod(b) );
    h2_hat =   h2_mod/h2_mod_modified;
    

%% detection using LLRT
    reconstruction_label3 = [h1_hat*ALP, h2_hat*ALP];
    [px] = optimalDetectorforN50(y.',8,reconstruction_label3.',1);
    px2 = sum(px(:,1:4),2);
    px3 = sum(px(:,5:8),2);
    
    for   k = 3:N/(fs)
        if sum(sum(log(px2(fs*(k-1)+1:fs*k,:)./px3(fs*(k-1)+1:fs*k,:))))>0
            c(k,:) = 0;
        else
            c(k,:) = 1;
        end
    end
    
%% convergence conditions    
    
     cha =   abs(h2_hat-h2_hat1)+abs(h1_hat-h1_hat1);  %The convergence conditions
 
    h1_hat1 = h1_hat;
    h2_hat1 = h2_hat;
       if cha< threshold
        break;  % If convergence, break the process
    end
end
end



function Px = optimalDetectorforN50(X,K,centroids, sigma)
%% The probability of data which follows Gussaian distribution
[N, D] = size(X);
j = sqrt(-1);
pMiu = centroids;
sigma_vec =  sigma*ones(D,1);
pSigma = diag(sigma_vec);
 Px = zeros(N, K);
        for k = 1:K
            Xshift = X-repmat(pMiu(k, :), N, 1);
            inv_pSigma = inv(pSigma);
            if rcond(pSigma) <= 3.079127e-18
                break;
            end
            tmp = real(sum((conj(Xshift)*inv_pSigma) .* Xshift, 2));
            coef = pi^(-D) * (det(pSigma)).^(-1);
            Px(:, k) = coef * exp(-tmp);
        end
end
