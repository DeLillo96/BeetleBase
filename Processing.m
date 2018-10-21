clear all; close all; clc;

tic

f0 = 3e9; %Central frequency equal to 3GHz
c = 3e8; %Light speed
dpsi = 2*pi/180; %Radar width in azimuth
dteta = 20*pi/180; %Radar width in elevation
eta = 0.5; %Antenna efficiency

lambda = c/f0; %Wavelength
r_max = 150e3; %Maximum distance from the Radar

r = [5e3:1e-1:r_max];

%A rectangular plane antenna is supposed:

L_a = lambda/dpsi; %Length in Azimuth
L_e = lambda/dteta; %Length in elevation
A = L_a*L_e; %Surface of the antenna

Gain = 4*pi*A*eta/(lambda^2); %Antenna gain

%We initially suppose a unitary transmitted power:

Pt = 1; %[W]

%Directive factor:

teta = 0;
f_dir = sinc(teta./dteta).^2;

%RCS of a common aircraft:

RCS = 50;

%Received power:

Pr = Pt./((4*pi*r.^2).^2)*Gain*A.*(f_dir.^2).*RCS;

figure();
plot(r/1e3,10*log10(Pr),'b','LineWidth',1.5),xlabel('r [km]'),grid on;
ylabel('P_{rx} [dB]'),title('Received power versus distance');

%The minimum received power is Pr = -185.8 dB. Thus the sensor must be
%dimensioned to be sensible to that amount of power received. However,
%the real received power is much greater since the transmitted power is
%much more than 1W: it will be shown later when we will 

PRI = 3*2*r_max/c; %Security factor equal to 3

%Band width choice:

B = 10e6; %10MHz or 100MHz
Toss = 30e-6; %Observation time

dt = 1/(2*B); %Scaling factor equal to 4
%Nt = 4*round(Toss/dt)+1; %Scegliamo un asse più largo per far vedere gli zeri
%t = [-(Nt-1)/2:(Nt-1)/2]*dt; 
t = [-PRI/2:dt:PRI/2];
Nt = length(t);
Nf = 2^ceil(log2(Nt)+2);

f = [-Nf/2:Nf/2-1]/Nf/dt;
df = f(2)-f(1);

%The transmitted signal is an up-chirp:

alfa = B/Toss;
g = exp(+1i*pi*alfa*(t).^2).*double(abs(t)<=Toss/2); %up-chirp

figure();
subplot(2,2,1),plot(t,real(g),'r'),grid on,title('Real part of g(t)');
xlabel('t [s]'),xlim([-1.5*Toss/2,1.5*Toss/2]);
subplot(2,2,3),plot(t,imag(g),'b'),grid on,title('Imaginary part of g(t)');
xlabel('t [s]'),xlim([-1.5*Toss/2,1.5*Toss/2]);
subplot(2,2,2),plot(t,abs(g),'r'),grid on,title('g(t) modulo');
xlabel('t [s]'),xlim([-1.5*Toss/2,1.5*Toss/2]);
subplot(2,2,4),plot(t,angle(g),'b'),grid on,title('g(t) phase');
xlabel('t [s]'),xlim([-1.5*Toss/2,1.5*Toss/2]);

%We want to guarantee a minimum SNR of 10 dB. In order to do so, we 
%choose that when an aircraft

teta_minPr = dteta/2;
f_dir_minPr = sinc((teta_minPr)./dteta).^2;
minPr = Pt./((4*pi*r_max.^2).^2)*Gain*A.*(f_dir_minPr.^2).*RCS;
Pr_rmax = Pt./((4*pi*r_max.^2).^2)*Gain*A.*(f_dir.^2).*RCS;

SNR = 10; %10dB = 10^(10/10) = 10

K = 1.38e-23; %Boltzmann constant
F = 3; %Noise factor
T0 = 290; %Temperature
N0 = K*T0*F; %Noise power spectral density

Eg = sum(abs(g).^2)*dt;
Pr_SNR = SNR*N0/Toss;
gamma = Pr_SNR/minPr; 
Pt = gamma*Pt; %This is the true transmitted power
G = fftshift(fft(g,Nf));
tau = abs(t(1));
G = G.*exp(1i*2*pi*f*tau)*dt; %Numerical delay compensation

figure();
subplot(2,1,1),plot(f/1e6,abs(G),'r'),grid on,title('G(f) modulo'),xlabel('f [MHz]');
subplot(2,1,2),plot(f/1e6,angle(G),'b'),grid on,title('G(f) phase'),xlabel('f [MHz]');

%Adapted filter:

h = conj(fliplr(g));

H = fftshift(fft(h,Nf));
H = H.*exp(1i*2*pi*f*tau)*dt;

%Signal after the adapted filter:

g1 = conv(g,h,'same')*dt; %Continuous convolution
G1 = G.*H;

figure();
subplot(2,1,1),plot(t,abs(g1),'r'),grid on,title('g''(t) modulo');
xlim([-1e-6,1e-6]),xlabel('t [s]');
subplot(2,1,2),plot(t,angle(g1),'b'),grid on,title('g''(t) phase');
xlim([-1e-6,1e-6]),xlabel('t [s]');

%Target simulation:

R = 150e3; %Aircraft distance [m]
V = 100; %Aircraft speed [m/s]

Pr = Pt./((4*pi*R.^2).^2)*Gain*A.*(f_dir.^2).*RCS;

%Received signal:

s_rx = 0;

for p = 1:length(R)
    tau = 2*R(p)/c;
    gp = double(abs(t-tau)<=Toss/2).*exp(1i*pi*alfa*(t-tau).^2);
    s_rx = s_rx + sqrt(Pr(p))*gp*exp(1i*4*pi/lambda*R(p));
end

%Noise simulation:

Pw = N0*B; %White noise with power spectral density equal to N0
w = sqrt(Pw/2)*(randn(1,Nt)+1i*randn(1,Nt));
%w = 0; %Remove % to not consider noise

s_rx = s_rx + w;

s_rc = conv2(s_rx,conj(fliplr(g)),'same')*dt;

% Peaks elevation:

Eg = sum(abs(g).^2)*dt; % = Toss
Peaks = sqrt(Pr)*Eg;

figure
r = t*c/2; %Convertion to spatial range
subplot(3,1,1),plot(t,abs(s_rx)),grid on;
xlabel('t [s]'),title('Received signal modulo')
subplot(3,1,2),plot(t,abs(s_rc)),grid on;
xlabel('t [s]'),title('Filtered signal modulo')
subplot(3,1,3),plot(r,abs(s_rc)),grid on,xlim([min(R)-1e2 max(R)+1e2]);
xlabel('d [m]'),title('Detail of the filtered signal modulo');

toc

%In order to evaluate the radial speed of the aircraft the Doppler effect
%must be taken into account. To that aim, the transmitted signal is 
%a signal composed by a down-chirp and an up-chirp in order not to have
%any ambiguity.

g = g + conj(g); %up-chirp + down-chirp

fd = 2/c*f0*V; %Doppler frequency

tau = 2*R/c;
g_rx = double(abs(t-tau)<=Toss/2).*...
        (exp(1i*pi*alfa*(t-tau).^2)+exp(-1i*pi*alfa*(t-tau).^2));
fase = 2*pi*f0*tau + 2*pi*fd*t;
s_rx = g_rx.*exp(-1i*fase);

df = 0.25/Toss;
f_max = min(0.5/dt,30*df);
f_dop_ax = (-f_max:df:f_max);
M = length(f_dop_ax);
s_rc_dop = zeros(M,length(t));

for m = 1:length(f_dop_ax)
    ss = s_rx.*exp(1i*2*pi*f_dop_ax(m)*t); %Frequencies shift
    s_rc_dop(m,:) = conv(ss,conj(fliplr(g)),'same')*dt;
end

r = c/2*t;
v = c/2*f_dop_ax/f0;
rho_r = c/2/B;
rho_v = c/2/f0/Toss;

figure
subplot(2,1,1), plot(t,abs(g_rx)), grid, ylim([0 max(abs(g_rx))*1.1])
xlabel('time [s]')
subplot(2,1,2), plot(f,abs(G)), grid, ylim([0 max(abs(G))*1.1])
xlabel('frequencies [Hz]')

figure
surf(r,v,abs(s_rc_dop)), 
shading interp, colormap('jet')
title(['Ambiguity function'])
xlabel('Range [m]'), ylabel('Speed [m/s]'),
xlim([0 r(end)])

figure
imagesc(r,v,abs(s_rc_dop)), axis xy, colormap('jet')
xlim([R-5*rho_r R+5*rho_r])
ylim([V-5*rho_v V+5*rho_v])
xlabel('Distance [m]')
ylabel('Speed [m/s]')
hold on, plot(R,V,'ko')
title('Ambiguity function - Maximum detail')