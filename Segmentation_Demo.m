%% Adaptive Segmentation Demo
% Create a synthetic signal and then segment it adaptively

%House Keeping Commands
clc
clear all
close all

%% Creates the Synthetic Signal and Visualizes its time varying structure
%Parameters
fs = 1024; % sampling frequency
t = 4096;% signal length (samples)
s1 = 256; %burst length
s2 = 640;%burst length
s3 = 128;%gap
s4 = 1664;%burst length
s5 = 1024; %burst length
s6 = 256;%gap
s7 = 128;%burst length

%-- Signal Components

% signal 1 
ed1=s1;
x1=[0:s1-1]/fs;
c1 = 2.*cos(2*pi*27*x1'); % 27Hz burst
xx1_ = c1.*(tukeywin(s1,0.75));
xx1 = [xx1_',zeros(1,(t-s1))];

%signal 2
ed2=(s1+s2);
x2 = [s1-128:ed2-1]/fs;
c2 = 1.05*cos(2*pi*23*x2'); % 22Hz burst
xx2_ = c2.*(tukeywin(s2+128,0.85));
xx2 = [zeros(1,s1-128),xx2_',zeros(1,t-ed2)];

% %%signal 3
 ed3=(ed2)+s3;
% x3 = [ed2-128:ed3-1]/fs;
% ns3= randn(length(x3),1)/20;
% xx3 = [zeros(1,ed2-128),ns3',zeros(1,t-ed3)];

%signal 4
ed4=(ed3+s4);
x4=[ed3-128:ed4-1]/fs;
c4 = .90*cos(2*pi*16*x4'); %16Hz burst
xx4_ = c4.*(tukeywin(s4+128,0.85));
xx4 = [zeros(1,ed3-128),xx4_',zeros(1,t-ed4)];

%signal 5
ed5=(ed4+s5);
x5=[ed4-128:ed5-1]/fs;
c5 = 1.5*cos(2*pi*25*x5'); % 25Hz burst
xx5_ = c5.*(tukeywin(s5+128,0.75));
xx5 = [zeros(1,ed4-128),xx5_',zeros(1,t-ed5)];

% %signal 6
 ed6=(ed5+s6);
% x6 = [ed5-128:ed6-1]/fs;
% ns6= randn(length(x6),1)/20;
% xx6 = [zeros(1,ed5-128),ns6',zeros(1,t-ed6)];

%signal 7
ed7=(ed6+s7);
x7=[ed6-128:ed7-1]/fs;
c7 = 1.75*cos(2*pi*21*x7'); %21Hz burst
xx7_ = c7.*(tukeywin(s7+128,0.85));
xx7 = [zeros(1,ed6-128),xx7_'];


%--Plot the Signal

x=[0:t-1]/fs; % signal length in seconds
%xx = xx1+xx2+xx3+xx4+xx5+xx6+xx7; %compose signal
xx = xx1+xx2+xx4+xx5+xx7; %compose signal

ns=randn(size(xx));
lfp = xx+ns/50; %add white noise to the entire signal

x=[0:length(lfp)-1]/fs;
set(figure(1),'pos',[281 552 595 199]); 
plot(x,lfp,'linewidth',1.2);
ylim([-3 3]);
xlabel('Time (s)'); ylabel('Amplitude');
set(gca,'FontSize',9,'Fontweight','bold','Linewidth',1)% changes axes numbers

%% Adaptive Segmentation of the Synthetic Signal

cns=generate_shifted_data(lfp',[-16:2:16]); %shifted datae
windw = 0.125; % min window length in seconds

% Construct Adaptive Segmentation
[ts,tsm,segind,et_seg,et_total]=lcpSeg_nondyadic(cns',windw,fs,1); % non-dyadic segmenting

set(figure(2),'pos',[281 552 595 199]);clf
pt=PlotTimeSeg_nondyadic(lfp',ts,fs);
xlabel('Time (s)'); ylabel('Amplitude');
set(gca,'FontSize',9,'Fontweight','bold','Linewidth',1)% changes axes numbers

%Time-frequency Map ----
[tf,t,f]=CreateTFMap2(lfp',ts,fs,2,64,2,1024);
set(figure(1),'pos',[282 259 646 198]); clf
imagesc(t,f,tf);
colorbar; axis xy;
ylim([ 0 40]); caxis([0 0.25]);

ylabel('Frequency (Hz)'); xlabel('Time (s)');
set(gca,'FontSize',9,'Fontweight','bold','Linewidth',1)% changes axes numbers
ts = [ts(1) diff(ts)];






