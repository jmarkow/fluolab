function [NEW_DATA,TIME]=fluolab_condition(DATA,TRIALS,FS,varargin)
%
%
%
%
% smooth data

blanking=[.2 0];
tau=.02;
newfs=100;
normalize=1;
detrend_win=.05;

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'FS'
			FS=varargin{i+1};
		case 'tau'
			tau=varargin{i+1};
		case 'newfs'
			newfs=varargin{i+1};
		case 'normalize'
			normalize=varargin{i+1};
		case 'detrend_win'
			detrend_win=varargin{i+1};
	end
end



decimate_f=round(FS/newfs);
cutoff=.8*(FS/2)/decimate_f;
cutoff=cutoff/(FS/2);

if ~isa(DATA,'double')
	DATA=double(DATA);
end

DATA=markolab_smooth(DATA,round(tau*FS));

[b,a]=ellip(3,.2,40,cutoff,'high');

NEW_DATA=filtfilt(b,a,DATA);
NEW_DATA=downsample(DATA,decimate_f);

[nsamples,ntrials]=size(NEW_DATA);

TIME=[1:nsamples]/newfs;
TIME=TIME(:);

NEW_DATA=fluolab_detrend(NEW_DATA,'fs',newfs,'win',detrend_win,'per',0);

if normalize==1
	tmp_max=max(NEW_DATA);
	tmp_min=min(NEW_DATA);

	tmp_min=repmat(tmp_min,[nsamples 1]);
	tmp_max=repmat(tmp_max,[nsamples 1]);

	NEW_DATA=(NEW_DATA-tmp_min)./(tmp_max-tmp_min);

elseif normalize==2

	NEW_DATA=zscore(NEW_DATA);

end


