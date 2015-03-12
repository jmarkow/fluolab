function C=fluolab_classify_trials(TTL,FS,varargin)
% check for noise/catch trials

if ~isa(TTL,'double')
	TTL=double(TTL);
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

noise_thresh=.6;
catch_thresh=.4;
noise_len=.06;
noise_pad=.05;
pulse_len=.02;
npulse=3;
interpulse_len=.02;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'noise_len'
			noise_len=varargin{i+1};
		case 'pulse_len'
			pulse_len=varargin{i+1};
		case 'npulse'
			npulse=varargin{i+1};
		case 'interpulse_len'
			interpulse_len=varargin{i+1};
	end
end

noise_len=round(noise_len*FS);
noise_pad=round(noise_pad*FS);
pulse_len=round(pulse_len*FS);
interpulse_len=round(interpulse_len*FS);

% construct detection filters

noise_filt=ones(noise_len,1);
noise_filt=[-ones(noise_pad,1);noise_filt;-ones(noise_pad,1)];
noise_filt=noise_filt/length(noise_filt);

catch_filt=[ones(pulse_len,1);-ones(interpulse_len,1)];
catch_filt=repmat(catch_filt,[npulse 1]);
catch_filt=catch_filt/length(catch_filt);

% get the noise trials

ttl_mat=double(TTL>1);
ttl_mat(ttl_mat==0)=-1;

noise_mat=filter(noise_filt,1,ttl_mat);
catch_mat=filter(catch_filt,1,ttl_mat);

[~,noise_trials]=find(noise_mat>noise_thresh);
noise_trials=unique(noise_trials);

[~,catch_trials]=find(catch_mat>catch_thresh);
catch_trials=unique(catch_trials);

C=zeros(1,size(TTL,2));
C(noise_trials)=1;
C(catch_trials)=2;
