function [precision, recall, fpr] = roc_stats(adj_true, adj_pred)
    p = size(adj_true,1);
    precision = sum((adj_true(:) ~= 0) & (adj_pred(:) ~= 0))/sum(adj_pred(:) ~= 0);
    recall = sum((adj_true(:) ~= 0) & (adj_pred(:) ~= 0))/sum(adj_true(:) ~= 0);
    fpr = sum((adj_true(:) == 0) & (adj_pred(:) ~= 0))/(sum(adj_true(:) == 0) - p);
end