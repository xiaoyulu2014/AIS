

%%% get index for worker
% if standard template:
n = str2num(getenv('SLURM_PROCID'));

% if job array template:
jobid = str2num(getenv('SLURM_ARRAY_JOB_ID'));
n = str2num(getenv('SLURM_ARRAY_TASK_ID'));

array_of_params = repelem(logspace(-4,0,20),10);

param=array_of_params(n+1); % param for this worker. Note need +1 since the ID's start on 0.
results = anis_1D(param); % all outputs to stdout

% if job array template:
results_file_string = ['/data/localhost/not-backed-up/xlu/results/anis/anis_1D' getenv('SLURM_ARRAY_JOB_ID') '_' getenv('SLURM_ARRAY_TASK_ID') '.mat']

save(results_file_string, 'results')
