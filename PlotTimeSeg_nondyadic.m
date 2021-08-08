function pt=PlotTimeSeg_nondyadic(eeg,ts,fs);
% pt=PlotTimeSeg(eeg,ts,fs);

maxeeg=max(max(abs(eeg)));
lts=length(ts);
maxts=max(ts(:));
t=(0:fs*maxts-1)/fs; %for regular zero to end time scale
y=ones(1,lts)*maxeeg*4;
Y=[y; -y];
T=repmat(ts,2,1); %for regular zero to end time scale
plot(t,eeg,'linewidth',1);
axis([0 maxts -1.5*maxeeg 1.5*maxeeg] );% for regular 0-end time scale
line(T,Y,'LineStyle','--','Color',[0 0 0]);
pt=1;
