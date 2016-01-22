function [RAW,REGRESS,TRIALS]=fluolab_fb_proc(DATA,AUDIO,TTL,varargin)
%  fiberdata
%
%
%
%
%

RAW=[];
REGRESS=[];
TRIALS=[];

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

channel=1;
daf_level=.3;
trial_cut=2;
normalize='m';
newfs=100;
dff=1;
tau=.1;
detrend_win=.3;
classify_trials='t';
detrend_method='p';
tau_regress=.05;
nmads=4;
detrend_sliding=0;
neg=0;
padding=[];
smooth_type='b';

for i=1:2:nparams
	switch lower(varargin{i})
		case 'daf_level'
			daf_level=varargin{i+1};
		case 'normalize'
			normalize=varargin{i+1};
		case 'classify_trials'
			classify_trials=varargin{i+1};
		case 'newfs'
			newfs=varargin{i+1};
		case 'dff'
			dff=varargin{i+1};
		case 'tau'
			tau=varargin{i+1};
		case 'detrend_win'
			detrend_win=varargin{i+1};
		case 'channel'
			channel=varargin{i+1};
		case 'detrend_method'
			detrend_method=varargin{i+1};
		case 'tau_regress'
			tau_regress=varargin{i+1};
		case 'nmads'
			nmads=varargin{i+1};
		case 'trial_cut'
			trial_cut=varargin{i+1};
		case 'neg'
			neg=varargin{i+1};
		case 'padding'
			padding=varargin{i+1};
    case 'smooth_type'
      smooth_type=varargin{i+1};
	end
end

if isempty(TTL) & isempty(AUDIO)
	classify_trials='';
end

if neg
	DATA.data=-DATA.data;
end

DATA=fluolab_datascrub(DATA,'channel',channel,'trial_cut',trial_cut,'nmads',nmads);

% where are the feedback trials?

[TRIALS,C]=fluolab_classify_trials(TTL,AUDIO,...
		'include_trials',DATA.trial_idx,'method',classify_trials,'daf_level',daf_level,...
		'padding',padding);

% code for change points, then window/ave

[RAW.mat,RAW.t]=fluolab_condition(DATA.data(:,:,channel),DATA.fs,DATA.t,'tau',tau,'detrend_win',detrend_win,...
	'newfs',newfs,'normalize',normalize,'dff',dff,'detrend_method',detrend_method,'smooth_type',smooth_type);
[nsamples,ntrials]=size(RAW.mat);

REGRESS.mat=markolab_deltacoef(RAW.mat',round(tau_regress*newfs),2)'; % approx 13 ms regression
REGRESS.t=RAW.t;

RAW.fs=newfs;
REGRESS.fs=newfs;

if nsamples==0
	warning('Data sample too short for analysis parameters...');
end

if ntrials<2
	warning('Not enough trials to continue...');
	return;
end

if any(isnan(RAW.mat(:)))
	warning('Found NaN in detrended data...');
end

% overall means

RAW=fluolab_mu(RAW,TRIALS);
REGRESS=fluolab_mu(REGRESS,TRIALS);
