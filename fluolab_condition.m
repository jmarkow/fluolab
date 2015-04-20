function [NEW_DATA,TIME]=fluolab_condition(DATA,FS,TIME,varargin)
% proc fluo data

tau=.02; % smoothing tau
newfs=100; % new sampling rate
normalize='m'; % normalize method
detrend_win=.05; % length of detrending window
dff=1; % dff?
detrend_method='p';

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'tau'
			tau=varargin{i+1};
		case 'newfs'
			newfs=varargin{i+1};
		case 'normalize'
			normalize=varargin{i+1};
		case 'detrend_win'
			detrend_win=varargin{i+1};
		case 'dff'
			dff=varargin{i+1};
		case 'detrend_method'
			detrend_method=varargin{i+1};
	end
end

% decimate data

decimate_f=round(FS/newfs);
cutoff=.8*(FS/2)/decimate_f;
cutoff=cutoff/(FS/2);

if ~isa(DATA,'double')
	DATA=double(DATA);
end

[nsamples,ntrials]=size(DATA);

[b,a]=ellip(3,.2,40,cutoff,'low');

NEW_DATA=filtfilt(b,a,DATA);
NEW_DATA=downsample(DATA,decimate_f);

[nsamples,ntrials]=size(NEW_DATA);
TIME=downsample(TIME,decimate_f);

tau_smps=round(tau*newfs);

if ~strcmp(lower(detrend_method(1)),'n')
	disp('Detrending...');
	NEW_DATA=fluolab_detrend(NEW_DATA,'fs',newfs,'win',detrend_win,'per',8,'dff',dff,'method',detrend_method);
end
NEW_DATA=markolab_smooth(NEW_DATA,tau_smps,'n','b');

NEW_DATA=NEW_DATA(tau_smps:end,:);
TIME=TIME(tau_smps:end);

[nsamples,ntrials]=size(NEW_DATA);


% exclude zero-padded portion at the beginning

% adjust time accordingly

switch lower(normalize(1))

	case 'm'

		tmp_max=max(NEW_DATA);
		tmp_min=min(NEW_DATA);

		tmp_min=repmat(tmp_min,[nsamples 1]);
		tmp_max=repmat(tmp_max,[nsamples 1]);

		NEW_DATA=(NEW_DATA-tmp_min)./(tmp_max-tmp_min);

	case 'z'

		NEW_DATA=zscore(NEW_DATA);

end


