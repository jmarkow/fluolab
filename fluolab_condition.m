function [NEW_DATA,TIME]=fluolab_condition(DATA,FS,varargin)
% proc fluo data

blanking=[.2 0]; % blanking interval
tau=.02; % smoothing tau
newfs=100; % new sampling rate
normalize='m'; % normalize method
detrend_win=.05; % length of detrending window
dff=1; % dff?

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'blanking'
			blanking=varargin{i+1};
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
	end
end

% decimate data

decimate_f=round(FS/newfs);
cutoff=.8*(FS/2)/decimate_f;
cutoff=cutoff/(FS/2);

if ~isa(DATA,'double')
	DATA=double(DATA);
end

DATA=markolab_smooth(DATA,round(tau*FS));
DATA=DATA(blanking(1)*FS:end-blanking(2)*FS,:);

[b,a]=ellip(3,.2,40,cutoff,'low');

NEW_DATA=filtfilt(b,a,DATA);
NEW_DATA=downsample(DATA,decimate_f);

[nsamples,ntrials]=size(NEW_DATA);

TIME=[1:nsamples]/newfs;
TIME=TIME(:);

NEW_DATA=fluolab_detrend(NEW_DATA,'fs',newfs,'win',detrend_win,'per',0,'dff',dff,'method','p');

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


