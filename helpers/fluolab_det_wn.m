function [ENERGY,BAND2_POWER,BAND1_POWER,T]=fluolab_det_wn(AUDIO,FS,varargin)
% check how flat the spectrum is each timepoint 
% regression?

if ~isa(AUDIO,'double')
	AUDIO=double(AUDIO);
end

if nargin<2
	disp('Setting FS to 30e3...');
	FS=30e3;
end

% the template cutoff could be defined by the 95th prctile of the abs(noise) magnitude

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs!');
end

len=.005; % window length (s)
band1=[1e3 5e3];
band2=[5e3 10e3];
overlap=0; % overlap (s)
ratio_thresh=2; % ratio song:nonsong
songpow_thresh=.8;
silence=0;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'song_band'
			song_band=varargin{i+1};
		case 'len'
			len=varargin{i+1};
		case 'overlap'
			overlap=varargin{i+1};
		case 'song_duration'
			song_duration=varargin{i+1};
		case 'ratio_thresh'
			ratio_thresh=varargin{i+1};
		case 'song_thresh'
			song_thresh=varargin{i+1};
		case 'pow_thresh'
			pow_thresh=varargin{i+1};
		case 'songpow_thresh'
			songpow_thresh=varargin{i+1};
		case 'silence'
			silence=varargin{i+1};

	end
end

len=round(len*FS);
overlap=round(overlap*FS);

[s,f,T]=spectrogram(AUDIO,len,overlap,[],FS);

% take the power and find our FS band

power=log(abs(s));

min_idx1=max(find(f<=band1(1)));
max_idx1=min(find(f>=band1(2)));

min_idx2=max(find(f<=band2(1)));
max_idx2=min(find(f>=band2(2)));

band1_power=mean(power(min_idx1:max_idx1,:));
band2_power=mean(power(min_idx2:max_idx2,:));

figure();plot(band1_power)
figure();plot(band2_power)
ENERGY=band2_power./band1_power;

% take the song/nonsong power ratio
