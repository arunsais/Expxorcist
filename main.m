function main(n, p, d, exNum, indir, outdir, base_measure_type, base_measure_args)
    rng(42);
    if ~isdeployed
        addpath(genpath([pwd '/algorithm/basis_expansion']));
        addpath(genpath([pwd '/stats']));
    end
    
    lambda = logspace(-1, 0.5, 18);
    
    precision_or = zeros(exNum, numel(lambda));
    recall_or = zeros(exNum, numel(lambda));
    fpr_or = zeros(exNum, numel(lambda));
    
    precision_and = zeros(exNum, numel(lambda));
    recall_and = zeros(exNum, numel(lambda));
    fpr_and = zeros(exNum, numel(lambda));
    
    train_nllk = cell(exNum, numel(lambda));
    test_nllk = cell(exNum, numel(lambda));
    
    params = cell(exNum, numel(lambda));
    
    infile = [indir '/data_' num2str(n) '_' num2str(p)];
    load(infile);
    adj = adj; % dummy line so that parfor works
    
    for k = 1 : exNum
        parfor l = 1 : numel(lambda)
            [params{k,l}, adj_pred_or, adj_pred_and, train_nllk{k,l}, test_nllk{k,l}] = npglm(xTrain{k}, xTest{k}, lambda(l), d, base_measure_type, base_measure_args);
            [precision_or(k,l), recall_or(k,l), fpr_or(k,l)] = roc_stats(adj, adj_pred_or);
            [precision_and(k,l), recall_and(k,l), fpr_and(k,l)] = roc_stats(adj, adj_pred_and);
            
            fprintf('%f, %f, %f, %f, %f, %f\n', recall_or(k,l), fpr_or(k,l), recall_and(k,l), fpr_and(k,l), train_nllk{k,l}(1), test_nllk{k,l}(1));
        end
    end
    
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end
    outfile = [outdir '/results_' num2str(n) '_' num2str(p) '_' num2str(d)];
    save(outfile, 'precision_or', 'recall_or', 'fpr_or', 'precision_and', 'recall_and', 'fpr_and', 'params', 'train_nllk', 'test_nllk');
    
    if ~isdeployed
        rmpath(genpath([pwd '/algorithm/basis_expansion']));
        rmpath(genpath([pwd '/stats']));
    end
end
