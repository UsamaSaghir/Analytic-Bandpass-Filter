%% Analytic filter demo
% Generetes a bandpass QAM signbal, designs Hilbert transformer, and completes analytic filtering for it.
% Normalized frequencies used, 1 corresponds to half sampling rate. 

% QAM test signal parameters:
sig='16QAM';    % Modulation order: QPSK/8PSK/16QAM/64QAM; defaiult '16QAM'
rolloff = 0.35; % Pulse-shaping rolloff; defualt 0.35
over = 16;      % Oversampling factor;   default 16
w_c = 0.4;      % Center frequency;      default 0.4

% Hilber filter parameters:
N_h = 30;       % FIR filter order
w_start = 0.1;  % Lower edge of passband
w_stop = 0.9;   % Higher edge of passband

fontti=12;       % Font size

if length(sig) == 4
    sig = [sig ' '];
end

if sig == 'QPSK '
	MM = 4;
	alphabet = exp(j*(0:MM-1)*2*pi/MM).*exp(j*pi/4);
elseif sig == '8PSK '
	MM = 8;
	alphabet = exp(j*(0:MM-1)*2*pi/MM).*exp(j*pi/8);
elseif sig == '16QAM'
	temp = [-3 -1 1 3];
	alphabet = repmat(temp,4,1) + j*repmat(temp,4,1).';
	alphabet = alphabet(:).';
elseif sig == '64QAM'
	temp = [-7 -5 -3 -1 1 3 5 7];
	alphabet = repmat(temp,8,1) + j*repmat(temp,8,1).';
	alphabet = alphabet(:).';
else
	error('Inconsistent input, aborting ...');
end

% index = randint(1,1000,[1 length(alphabet)]);			% Random QAM symbol sequence, signal
index = randi(length(alphabet),1,1000);			% Random QAM symbol sequence, signal
sym = alphabet(index);

pulse = rcosine(1,over,'normal',rolloff);	% basic raised-cosine pulse-shape
[trash,pos] = max(pulse);				   % filter delay 'pos'

sig = kron(sym,[1 zeros(1,over-1)]);	% oversampled signal, see help kron
sig = filter(pulse,1,sig);				   % signal after pulse shaping
sig = sig(pos:length(sig));				% discard transient

r = 2*real(sig.*exp(j*w_c*pi*(0:length(sig)-1))); 	% Bandpass signal

W = linspace(-pi,pi,1024);			% Discrete-time frequency axis
R = freqz(r(1000:2023),1,W); 		% Spectrum of the passband signal
R = R./max(abs(R));			      % Normalize the spectrum for plotting

% Hilbert filter 
h = remez(N_h,[w_start w_stop],[1 1],'Hilbert');	      % Design of a wideband Hilbert filter of order N_h	
pos = N_h/2 + 1;					                           

f2 = 0.5*([zeros(1,N_h/2) 1 zeros(1,N_h/2)] + j*h);	
F2 = freqz(f2,1,W);

y = filter(h,1,r);				% Signal filtering
y = y(pos:length(y));			% discard transient 

y = r(1:length(y))+j*y;

Y = freqz(y(1000:2023),1,W);	% Spectrum of the signal at this stage
Y = Y./max(abs(Y));				% Normalize the spectrum for ploting

% figure; 
% hFig1 = figure('visible', 'off');
figure

subplot(211);                           % Plot again the bandpass signal entering the receiver
plot(W/pi,abs(R));grid;
hold on;
plot(W/pi,abs(F2),'r','LineWidth',2);
set(gca,'FontSize',fontti);
title('Signal spectrum and analytic filter frequency response');
ylabel('Amplitude','FontSize',fontti);
axis([-1 1 0 1.1]);

hold off;

subplot(212);
plot(W/pi,abs(Y));grid;
set(gca,'FontSize',fontti);
title('Signal after analytic filtering');
xlabel('Frequency \omega / \pi','FontSize',fontti);
ylabel('Amplitude','FontSize',fontti);
axis([-1 1 0 1.1]);

%Render jpegs and write to files.
drawnow;
% s.GraphFileName1 = sprintf('%sas1.jpeg', mlid); % Change the name to your image name
% wsprintjpeg(hFig1, s.GraphFileName1);
% s.GraphFileName1 = sprintf('/invocom/p3-1/%sas1.jpeg',mlid); % Change the path to your path
% 
% close(hFig1)
% 
% % Put name of graphic file into HTML template file.
% templatefile = which('as.htm');
% rs = htmlrep(s, templatefile);	 
