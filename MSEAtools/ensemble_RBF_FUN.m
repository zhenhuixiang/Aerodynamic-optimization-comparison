% tool function ensemble_RBF_FUN 
function y = ensemble_RBF_FUN(x, RBF_FUN_Arr)
    fit = [];
    model_num = length(RBF_FUN_Arr);
    for i = 1:model_num
        RBF_FUN = RBF_FUN_Arr{i};
        fit(i,:) = RBF_FUN(x);
    end
    y = mean(fit,1);
end