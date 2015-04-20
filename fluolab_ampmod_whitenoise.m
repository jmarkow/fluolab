function AMPMOD_WN=fluolab_ampmod_whitenoise(SONG,varargin)
%
%
%
%
%

filtering=300;
fs=24414.0625;
norm_amp=1;
padding=[.2 .2];
nfft=[];

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'filtering'
			filtering=varargin{i+1};
		case 'norm_amp'
			norm_amp=varargin{i+1};
		case 'padding'
			padding=varargin{i+1};
		case 'nfft'
			nfft=varargin{i+1};
	end
end

padding=round(padding*fs);
SONG=SONG(padding(1):end-padding(2),:);

if ~isa(SONG,'double')
	SONG=double(SONG);
end

if ~isempty(filtering)
	[b,a]=ellip(5,.2,40,[filtering]/(fs/2),'high');
	SONG=filtfilt(b,a,SONG);
end

if isvector(SONG)
	SONG=SONG(:);
end

[nsamples,ntrials]=size(SONG);

if isempty(nfft)
	nfft=nsamples;
end

% assume song is samples x trials

% get random phases
%

scr_theta=angle(fft(rand(nfft,1)));

sig_fft=fft(SONG,nfft);

if ntrials>1
	sig_fft=mean(sig_fft,2);
end

sig_amp=abs(sig_fft);
AMPMOD_WN=real(ifft(sig_amp.*exp(1j.*scr_theta)));

if norm_amp
	AMPMOD_WN=AMPMOD_WN./max(abs(AMPMOD_WN));
end
