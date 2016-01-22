function DATA=fluolab_mu(DATA,TRIALS,varargin)
%
%
%
%

nboots=1e3;

nparams=length(varargin);
if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'nboots'
			nboots=varargin{i+1};
  end
end

trial_types=fieldnames(TRIALS.fluo_include);
ntypes=length(trial_types);

for i=1:ntypes

	if length(TRIALS.fluo_include.(trial_types{i}))<3

		DATA.ci.(trial_types{i})=[];
		DATA.mu.(trial_types{i})=[];

		continue;

	end

	if strcmp(trial_types{i},'include'), continue; end
	if strcmp(trial_types{i},'idx'), continue; end

	DATA.ci.(trial_types{i})=bootci(nboots,{@mean,DATA.mat(:,TRIALS.fluo_include.(trial_types{i}))'},'type','cper');
	DATA.mu.(trial_types{i})=mean(DATA.mat(:,TRIALS.fluo_include.(trial_types{i}))');

end
