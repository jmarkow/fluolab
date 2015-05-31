function [TRIALS,C]=fluolab_classify_trials(TTL,AUDIO,varargin)
% check for noise/catch trials



if nargin<2
	AUDIO=[];
end

if isstruct(AUDIO) & ~isa(AUDIO.data,'double')
	AUDIO.data=double(AUDIO.data);
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
blanking=[.15 .15];
method='t'; % [t]tl or [s]ound
daf_level=.4;
smoothing=.06;
daf_cutoff=500;

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
		case 'noise_thresh'
			noise_thresh=varargin{i+1};
		case 'catch_thresh'
			catch_tresh=varargin{i+1};
		case 'blanking'
			blanking=varargin{i+1};
		case 'method'
			method=varargin{i+1};
		case 'daf_level'
			daf_level=varargin{i+1};
		case 'smoothing'
			smoothing=varargin{i+1};
	end
end


switch lower(method(1))

	case 't'

		if ~isa(TTL.data,'double')
			TTL.data=double(TTL.data);
		end

		%C=fluolab_classify_trials(TTL.data(blanking_idx,:),TTL.fs);

		[nsamples,ntrials]=size(TTL.data);

		noise_len=round(noise_len*TTL.fs);
		noise_pad=round(noise_pad*TTL.fs);
		pulse_len=round(pulse_len*TTL.fs);
		interpulse_len=round(interpulse_len*TTL.fs);

		% construct detection filters

		blanking_idx=[];

		if blanking(1)==0
			blanking_idx(1)=1;
		else
			blanking_idx(1)=round(blanking(1)*TTL.fs);
		end

		if blanking(2)==0
			blanking_idx(2)=nsamples;
		else
			blanking_idx(2)=nsamples-round(blanking(2)*TTL.fs);
		end

		blanking_idx=blanking_idx(1):blanking_idx(2);

		noise_filt=ones(noise_len,1);
		noise_filt=[-ones(noise_pad,1);noise_filt;-ones(noise_pad,1)];
		noise_filt=noise_filt/length(noise_filt);

		catch_filt=[ones(pulse_len,1);-ones(interpulse_len,1)];
		catch_filt=repmat(catch_filt,[npulse 1]);
		catch_filt=catch_filt/length(catch_filt);

		% get the noise trials

		ttl_mat=double(TTL.data(blanking_idx,:)>1);
		ttl_mat(ttl_mat==0)=-1;

		noise_mat=filter(noise_filt,1,ttl_mat);
		catch_mat=filter(catch_filt,1,ttl_mat);

		[~,noise_trials]=find(noise_mat>noise_thresh);
		noise_trials=unique(noise_trials);

		[~,catch_trials]=find(catch_mat>catch_thresh);
		catch_trials=unique(catch_trials);

		C=zeros(1,size(TTL.data,2));
		C(catch_trials)=2;
		C(noise_trials)=1; % noise_trial catch_trial overlap overwritten as noise trial

		TRIALS.all.catch=find(C==2);
		TRIALS.all.daf=find(C==1);
		TRIALS.all.other=find(C==0);
		TRIALS.all.catch_other=find(C==0|C==2);

	case 's'

		if ~isa(AUDIO.data,'double')
			AUDIO.data=double(AUDIO.data);
		end

		[nsamples,ntrials]=size(AUDIO.data);
	
		blanking_idx=[];

		if blanking(1)==0
			blanking_idx(1)=1;
		else
			blanking_idx(1)=round(blanking(1)*AUDIO.fs);
		end

		if blanking(2)==0
			blanking_idx(2)=nsamples;
		else
			blanking_idx(2)=nsamples-round(blanking(2)*AUDIO.fs);
		end

		blanking_idx=blanking_idx(1):blanking_idx(2);

		[b,a]=ellip(4,.2,40,[daf_cutoff]/(AUDIO.fs/2),'high');
		daf_mat=filtfilt(b,a,AUDIO.data(blanking_idx,:)).^2;
		smooth_smps=round(smoothing*AUDIO.fs);
		daf_mat=filter(ones(smooth_smps,1)/smooth_smps,1,daf_mat)>daf_level;
		
		[~,daf_trials]=find(daf_mat);
		daf_idx=unique(daf_trials);
		first_daf=1;

		daf_trials=zeros(ntrials,1);
		daf_trials(daf_idx)=1;

		C=daf_trials;

		% for now power detection 10e3-12e3 works fine

		% catch trials

		TRIALS.all.catch=[];
		TRIALS.all.daf=find(daf_trials==1);
		TRIALS.all.other=find(daf_trials==0);
		TRIALS.all.catch_other=[];

	otherwise

end


