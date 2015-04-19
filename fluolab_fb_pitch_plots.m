function PITCHFIGS=fluolab_fb_plots(PITCH,varargin)


nparams=length(varargin);
visible='on';
datenums=[];
smooth_trials=100;
bin_res=5;
hist_colors='winter';
hist_order=1e2;
ylim_order=1e2;
pitch_threshold=[];
pitch_condition='';

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'datenums'
			datenums=varargin{i+1};	
		case 'visible'
			visible=varargin{i+1};
		case 'smooth_trial'
			smooth_trials=varargin{i+1};
		case 'hist_order'
			hist_order=varargin{i+1};
		case 'pitch_threshold'
			pitch_threshold=varargin{i+1};
		case 'pitch_condition'
			pitch_condition=varargin{i+1};
	end
end

if ~isempty(pitch_condition) & (strcmp(lower(pitch_condition(1)),'l') | strcmp(lower(pitch_condition(1)),'g'))
	iscondition=1;
else
	iscondition=0;
end

isdate=~isempty(datenums);
multdate=iscell(datenums);
multpitch=iscell(PITCH);
isthresh=~isempty(pitch_threshold);

hist_colors=colormap([hist_colors '(' num2str(length(PITCH)) ')']);

% plot pitch from today

if isdate & ~multpitch

	% timecourse
	
	trial_mu=mean(PITCH);
	
	ylimits=prctile(trial_mu,[5 95])
	ylimits=ylimits/ylim_order;
	ylimits(1)=floor(ylimits(1))*ylim_order;
	ylimits(2)=ceil(ylimits(2))*ylim_order;
	
	[nsteps,ntrials]=size(PITCH);

	PITCHFIGS.pitch_timecourse=figure('paperpositionmode','auto','visible',visible);

	if isthresh & iscondition

		if strcmp(lower(pitch_condition(1)),'g')

			ydata(1,:)=repmat(pitch_threshold,[1 length(datenums)]);
			ydata(2,:)=repmat(ylimits(2),[1 length(datenums)]);
		
			facecolor=[.5 1 1];

		else
			ydata(1,:)=repmat(ylimits(1),[1 length(datenums)]);
			ydata(2,:)=repmat(pitch_threshold,[1 length(datenums)]);

			facecolor=[1 1 .5];

		end

		markolab_shadeplot(datenums,ydata,facecolor,'none');
		hold on;

	end

	plot(datenums,trial_mu,'k.','markersize',15);
	
	% smoothed timecourse

	smooth_pitch=smooth(datenums,trial_mu,min(smooth_trials,ntrials));	

	hold on;
	plot(datenums,smooth_pitch,'b-');
	datetick('x','HH:MM (PM)');
	box off;
	ylabel('Pitch (Hz)');


	ylim([ylimits]);
	xlim([datenums(1) datenums(end)])
	set(gca,'layer','top');

	% show summary stats

end

if multpitch & multdate & length(multdate)==length(multpitch)

	% cycle through cell arrays, plot entire time-course, etc. etc.
	
	PITCHFIGS.multpitch_timecourse=figure('paperpositionmode','auto','visible',visible);
	
	for i=1:length(PITCH)

		trial_mu=mean(PITCH{i});
		[nsteps,ntrials]=size(PITCH{i});

		days=daysdif(datenums{1}(1),datenums{i});

		plot(days,trial_mu,'k.');
		hold on;

		% smoothed timecourse

		smooth_pitch=smooth(days,trial_mu,min(smooth_trials,ntrials));
		plot(days,smooth_pitch,'b-');

	end

	box off;
	ylabel('Pitch (Hz)');

	tmp=mean(cat(2,PITCH{:}))

	ylimits=prctile(tmp,[5 95]);
	ylimits=ylimits/ylim_order;
	ylimits(1)=floor(ylimits(1))*ylim_order;
	ylimits(2)=ceil(ylimits(2))*ylim_order;
	ylim([ylimits]);


end


if multpitch

	% histogram everything

	tmp=mean(cat(2,PITCH{:}));
	bin_edges=prctile(tmp,[5 95]);
	bin_edges=(bin_edges/hist_order);
	bin_edges(1)=floor(bin_edges(1))*hist_order;
	bin_edges(2)=ceil(bin_edges(2))*hist_order;

	bins=bin_edges(1):bin_res:bin_edges(2);

	PITCHFIGS.multpitch_hist=figure('paperpositionmode','auto','visible',visible);

	for i=1:length(PITCH)

		n=histc(mean(PITCH{i}),bins);
		n=n./sum(n);

		markolab_stairplot(n,bins,'color',hist_colors(i,:));
		hold on;

	end

end


if ~multpitch
	PITCHFIGS.pitch_hist=figure('paperpositionmode','auto','visible',visible);

	tmp=mean(PITCH);
	bin_edges=prctile(tmp,[5 95]);
	bin_edges=(bin_edges/hist_order);
	bin_edges(1)=floor(bin_edges(1))*hist_order;
	bin_edges(2)=ceil(bin_edges(2))*hist_order;

	bins=bin_edges(1):bin_res:bin_edges(2);
	n=histc(tmp,bins);
	n=n./sum(n);
	markolab_stairplot(n,bins,'color',hist_colors(1,:));
	xlabel('Pitch (Hz)');
	ylabel('P');
	box off;

end
