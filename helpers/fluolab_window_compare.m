function [win_data]=fluoflab_sliding_window(DATA,TTL,TRIALS,varargin)
% plot window centered on ttl pulse
%
%
%
%
colors='jet';
cbar_dist=.1;
cbar_width = .025; 
clim_order=1e3;
datenums=[];
time_dist=.005;
time_width=.075;
escape_dist=.005;
escape_width=.075;
escape_idx=[];
time_win=[.1 .2];
padding=[.2 .2];
thresh=.5;

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'cbar_dist'
			cbar_dist=varargin{i+1};
		case 'colors'
			colors=varargin{i+1};
		case 'cbar_width'
			cbar_width=varargin{i+1};
		case 'datenums'
			datenums=varargin{i+1};
		case 'label'
			label=varargin{i+1};
		case 'escape_idx'
			escape_idx=varargin{i+1};
		case 'padding'
			padding=varargin{i+1};
        case 'time_win'
            time_win=varargin{i+1};
	end
end

[nsamples,ntrials]=size(DATA.mat);
TTL.data=TTL.data(:,TRIALS.all.fluo_include);

% how to estimate mean/ci

mufun=@(x) mean(x,2);
varfun=@(x) std(x,[],2);

% get onsets

[nsamples_ttl]=size(TTL.data,1);

idx=(1:nsamples_ttl-1);

% pos_cross

pos_cross=zeros(1,ntrials);
data_fs=DATA.fs;

win_len=round(time_win.*DATA.fs);
tmp=zeros(sum(win_len)+1,ntrials);

padding_smps=round(padding.*DATA.fs);
todel=[];

for i=1:ntrials
	
	cur_ttl=TTL.data(:,i);
	i
	hit=TTL.t(find(cur_ttl(idx)<thresh&cur_ttl(idx+1)>=thresh));
	hit(hit<padding(1)|hit>TTL.t(end)-padding(2))=[];

	% map the hit to the nearest timepoint in the singing data

	if isempty(hit)
		todel=[todel i];
		continue;
	end

	hit=hit(1)
	[~,loc]=min(abs(DATA.t-hit));
	
	pos_cross(i)=loc;
    
	if loc-win_len(1)<=0 | loc+win_len(2)>nsamples
		todel=[todel i];
		continue;
    end
    
	tmp(:,i)=DATA.mat(loc-win_len(1):loc+win_len(2),i);

end

todel
% distribute tmp

win_data.daf=tmp(:,setdiff(TRIALS.fluo_include.daf,todel));
win_data.catch=tmp(:,setdiff(TRIALS.fluo_include.catch,todel));
win_data.t=[-win_len(1):win_len(2)]'/DATA.fs;


