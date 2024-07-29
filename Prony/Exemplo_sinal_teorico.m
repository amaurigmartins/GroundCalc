clear all
close all
clc

%      Exemplo de sinal teórico para aplicação da análise Prony:
%
% o sinal é composto por uma componente de frequência fundamental, uma com-
% ponente CC de decaimento exponencial e duas componentes de frequência
% sub-síncronas. O sinal é adicionado de ruído com distribuição uniforma
% com média 0 e desvio padrão de 1% da freq. fundamental e 2 - 20
% harmônicos com 5% de THD

N = 64;
f = 60;
fs = N*f;
w1 = 2*pi*f;

t = 0:1/fs:4/f-1/fs;

A1      = 3;
B_DDOC  = 0.23;
B_SSFC1 = 4.4;
B_SSFC2 = 8;

teta1      = -47*pi/180;
teta_SSFC1 = 142*pi/180;
teta_SSFC2 = -17*pi/180;

f_SSFC1 = 30;
f_SSFC2 = 42;

tau_0 = -0.04;
tau_SSFC1 = -0.024;
tau_SSFC2 = -0.02;

std_f1 = 0.01*60; % Desvio padrão de 1% de 60 Hz

comp_fund  = A1*cos(w1*t + teta1);

comp_DDOC  = B_DDOC*exp(t/tau_0);

comp_SSFC1 = B_SSFC1.*exp(t/tau_SSFC1).*cos(2*pi*f_SSFC1*t + teta_SSFC1);

comp_SSFC2 = B_SSFC2.*exp(t/tau_SSFC2).*cos(2*pi*f_SSFC2*t + teta_SSFC2);

noise = sqrt(std_f1)*randn(size(t));

i = comp_fund + comp_DDOC + comp_SSFC1 + comp_SSFC2 + noise;

% Harmônicas de 2ª até 20ª componentes com THD de 5%
for k = 2:20
    i = i + 0.034*cos(k*w1*t);
end

figure(1)
plot(t,i)
title('Sinal teórico')
xlabel('Tempo')
ylabel('Amplitude')


% Passando o filtro de média e descartando já o primeiro ciclo após a falta

i_med = zeros(1,length(t));

for n = N+1:length(t)
    for k = 1:N
        i_med(n) = i_med(n) + i(n-k);
    end
    i_med(n) = (1/N)*i_med(n) + (i(n)+i(n-N+1))/(2*N);
end

figure(2)
plot(t(N:length(t)),i(N:length(t)),t(N:length(t)),i_med(N:length(t)))
title('Sinal obtido do filtro de média')
xlabel('Tempo')
ylabel('Amplitude')
legend('Teórico','Média')


ini = N+1;
fim = length(t);
tx  = t(ini:fim)-t(ini);

i_medx = i_med(ini:fim);

figure(3)
plot(t,i_med,tx+(N)/fs,i_medx),legend('i_med','i_medx')


% [xData, yData] = prepareCurveData( tx, i_medx );
% 
% % Set up fittype and options.
% ft = fittype( 'fourier6' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.Normalize = 'on';
% opts.StartPoint = [0 0 0 0 0 0 0 0 0 0 0 0 0 1.21864237941788];
% 
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% 
% % Plot fit with data.
% figure( 'Name', 'Sinal teórico' );
% h = plot( fitresult, xData, yData );
% legend( h, 'i_medx vs. tx', 'Sinal teórico', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'tx', 'Interpreter', 'none' );
% ylabel( 'i_medx', 'Interpreter', 'none' );
% 
% figure(4)
% plot(t,i_med,tx+(N)/fs,i_medx,xData+(N)/fs,yData),legend('i_med','i_medx','Fit')


% Definição do sinal transitório
% i_medx_fit(1:N) = 0;
% i_medx_fit(N+1:length(t)) = yData;
% 
% figure(5)
% plot(t,i_med,t,i_medx_fit),legend('i_med','Sinal fit')



% Aplicação prony

test_data = i_medx;
test_time = tx;
N_order = 90; % ordem
SUB_N = 10;
criteria_val = 1;

test_t_increment=tx(2)-tx(1);
fs=1/test_t_increment; % Sampling frequency

[inumz,idenz]=prony(test_data,N_order,N_order);

% Compute and plot the impulse response of the Prony approximating system
[iapp,tapp]=impz(inumz,idenz,length(test_time),fs);
% [r,p,k]=residuez(inumz,idenz);
% a_list=r(:);
% spoles=log(p(:))*fs;
% tau_list=1./real(spoles);
% omega_list=imag(spoles);
% 
% n_size=size(test_time,1); 
% n=0:1:n_size;
% for ct=1:N_order
%     energy(ct)=abs(a_list(ct)^2)*(sum(abs(p(ct).^n).^2));
% end
% if (criteria_val==1)
%     [Mag, ISort]=sort(abs(a_list));
% else
%     [En,ISort]=sort(energy);
% end
% ISort=ISort(end:-1:1);
% FULL_IND=[ISort(:)]' ;
% ai=zeros(size(test_time));
% SUB_IND=FULL_IND(1:SUB_N);
% for i=SUB_IND
%    ai=ai+a_list(i)*exp(spoles(i)*tapp);
% end
% ai=real(ai);
% % Remove bug, PA should be same if modes are equal to model order
% if(SUB_N==N_order)
%     ai(1)=iapp(1);
% end
% % Update the SUB_N for prony fit when no. of modes is not specified
% if(SUB_N==0)
%     SUB_IND=FULL_IND;
% end
% 
figure(6)
plot(tapp,iapp,tx,i_medx),legend('Prony','I_med')



%% VE SE ISSO PRESTA

x=test_data.';
p=192/2;
Ts= test_time(end)-test_time(end-1);
[Amp,alfa, freq, theta] = polynomial_method(x,p,Ts,'tls')

[freq, ix] = sort(freq);
Amp = Amp(ix);
alfa = alfa(ix);
theta = theta(ix);

% maximum number of components (positive frequencies)
MAX_COMP = length(find(freq>=0));

ii = find(freq<0,1,'last');


Ncomp = MAX_COMP;


% get the higher positive frequency (fcut) and filter from its
% negative value (-fcut) to the positive one (fcut)

fcut = freq(ii+Ncomp); % fcut from 0Hz to selected component number
ix = find(freq>=-fcut & freq<=fcut);

comp = length(ix); % number of components used to recover the signal

% filtering step
freq_ = freq(ix);
Amp_ = Amp(ix);
alfa_ = alfa(ix);
theta_ = theta(ix);

% recovering temporal filtered signal

N=length(x);

n = 0:N-1;
n = repmat(n, [comp 1]);
k = exp(repmat(alfa_, [1 N]).*n*Ts);

k(isinf(k))=realmax*sign(k(isinf(k)));

Y = transpose(repmat(Amp_,[1 N]).*k)*...
cos(2*pi*Ts*repmat(freq_, [1 N]).*n + repmat(theta_,[1 N]));

yy = diag(Y);

figure(7);clf
plot(tapp,yy,'-o');hold all
plot(tx,i_medx)
legend(sprintf('Prony order %d ',p), 'Imed')