function [DATA]=fluolab_fb_proc_window(DATA,TTL,varargin)
%  fiberdata
%
%
%
%
%

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

blanking=[.05 .05];
channel=1;
ttl_level=1;
trial_cut=2;
normalize=1;
win_size=[.05 .2];

for i=1:2:nparams
	switch lower(varargin{i})
		case 'blanking'
			blanking=varargin{i+1};
		case 'daf_level'
			daf_level=varargin{i+1};
		case 'normalize'
			normalize=varargin{i+1};
	end
end

% where are the feedback trials?

[nsamples,ntrials]=size(TTL.data);
blanking_idx=[round(blanking(1)*TTL.fs):nsamples-round(blanking(2)*TTL.fs)];

daf_mat=(TTL.data(blanking_idx,:)>ttl_level);
[smps,trials]=find(daf_mat);
uniq_trials=unique(trials)
ntrials=length(uniq_trials)

win_smps=round(win_size*DATA.fs);
win_vec=[-win_smps(1):win_smps(2)];
win_len=length(win_vec)

DATA.win.mat=zeros(win_len,trials);

ttl_nsamples=size(TTL.data,1);
data_nsamples=size(DATA.mat,1);

ttl_t=[1:ttl_nsamples]/TTL.fs;
data_t=[1:data_nsamples]/DATA.fs;

to_del=[];

% need to separate daf from catch!

for i=1:ntrials

	% get the rise time
	
	daf_samples=smps(trials==uniq_trials(i));

	% first sample
	
	first_t=ttl_t(min(daf_samples));

	% convert to fluoresence timebase
	
	[~,idx]=min(abs(first_t-data_t));

	% get the window if we can


	startpoint=idx-win_smps(1);
	stoppoint=idx+win_smps(2);

	if startpoint<1 | stoppoint>data_nsamples
		to_del=[to_del i];
		continue;
	end

	DATA.win.mat(:,uniq_trials(i))=DATA.mat(idx-win_smps(1):idx+win_smps(2),uniq_trials(i));
end

