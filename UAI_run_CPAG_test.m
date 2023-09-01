function UAI_run_CPAG_test(sessionID, NN, d, p, nLoop)
% UAI batch run of nLoop cg-to-cpag tests over different graph sizes (in NN
% array), with density d and cycle prob. p
% output is saved as mat-file 'datafile', processed by 'analRes_UAIrun'

% 1 - configure runtest parameters
global DBG
DBG = 0;
if ~isempty(DBG), DEBUG = DBG; else DEBUG = 0; end;

fprintf('\nRunning session %s, [1:%d] runs over graphs size %d-%d for density %d\n', ...
          sessionID, nLoop, NN(1),NN(end), d);

datafile = sprintf('Test_UAI_CPAG_%s_pC%d.mat',sessionID, nLoop);
fprintf('datafile = %s \n',datafile);
           
% ===============================================
% initialize
ACC   = [];

% 2 - create structure mapping (initialize once)

fprintf('\n ================== \n RUNNING %d loops:\n',nLoop);

nN = length(NN);
% =========================================
for ni = 1:nN
  % loop over graphs size N
  N  = NN(ni);
    
  % repeat nLOOP times
  for entry = 1:nLoop
    % get entry in ACC array
    idx = (ni-1)*max(nLoop) + entry;
    % initialize rng per graph (for reproducability)
    state = idx*20+101;
    rng(state);

    % generate model-set
    time = clock;
    fprintf('[%d:%d:%2.1f] === RUNNING %d nodes, loop %d:\n',time(4:6),N,entry);

    % run algorithms
    % create graph
    [G,S,A,M] = mk_random_cg(N,d,p);
    % get nr |SCC|>= 2, max
    SCC   = sum(S,1);
    maxS  = max(SCC);
    nrSCC = size(S,2);
    nrCyc = sum(SCC > 1);

    % run original CPAG-from-Graph
    t0 = tic;
    [P,vT,nCounts] = cpag_from_cg_org(G);
    nCounts(3,end) = toc(t0);       % track total time spent
    nrMark_1   = sum(sum(P > 0));
    nrOrient_1 = sum(sum(P == 1 | P == 2));
    nrvT_1     = size(vT,1);
    
    % run new Graph-to-CPAG
    t0 = tic;
    [P2,vT2,uT2,M,nCounts2] = cg_to_cpag_new(G);
    nCounts2(3,end) = toc(t0);      % track total time spent
    nrMark_2   = sum(sum(P2 > 0));
    nrOrient_2 = sum(sum(P2 == 1 | P2 == 2));
    nrvT_2     = size(vT2,1);

    % store result
    ACC.Info(idx, : ) = [N,d,entry,maxS,nrSCC,nrCyc,(nrMark_1 - nrMark_2), ...
                         (nrOrient_1 - nrOrient_2), (nrvT_1 - nrvT_2)];
    ACC.R1(idx , :) = [N,d,entry,nCounts(:)'];       % so end = 3+21 = 24
    ACC.R2(idx , :) = [N,d,entry,nCounts2(:)'];

    % save result (yes each time overwrite file ...
  end;  % for entry
  save(datafile, 'ACC');
end;  % for NN

disp('Finished');


  