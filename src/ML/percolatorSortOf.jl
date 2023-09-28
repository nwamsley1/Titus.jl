function getQvalues!(PSMs::DataFrame, probs::Vector{Union{Missing, Float64}}, labels::Vector{Union{Missing, Bool}})
    #Could bootstratp to get more reliable values. 
    q_values = zeros(Float64, (length(probs),))
    order = reverse(sortperm(probs)) #Sort class probabilities
    targets = 0
    decoys = 0
    for i in order
        if labels[i] == true
            decoys += 1
        else
            targets += 1
        end
        q_values[i] = decoys/(targets + decoys)
    end
    PSMs[:,:q_value] = q_values;
end

getQvalues!(PSMs::DataFrame, probs::Vector{Float64}, labels::Vector{Bool}) = getQvalues!(PSMs, allowmissing(probs), allowmissing(labels))

function rankPSMs!(PSMs::DataFrame, features::Vector{Symbol}; n_folds::Int = 3, colsample_bytree::Float64 = 0.5, num_round::Int = 25, eta::Float64 = 0.15, min_child_weight::Int = 1, subsample::Float64 = 0.5, gamma::Int = 0, max_depth::Int = 10, max_train_size::Int = 1e6, print_importance::Bool = false)
   
    #[:hyperscore,:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all, :spectral_contrast_matched,:RT_error,:scribe_score,:y_ladder,:b_ladder,:RT,:diff_hyper,:median_ions,:n_obs,:diff_scribe,:charge,:city_block,:matched_ratio,:weight,:intensity,:count,:SN]
    X = Matrix(PSMs[:,features])
    X_labels = PSMs[:, :decoy]

    #Using a random selection of rows means that a pair strongly correlated rows (adjacent scans) could end up
    #split between the training and testing fold, potentially causing optimistic error estimates. Perhaps
    #this sampling scheme should be changed in the future?
    permutation = randperm(size(PSMs)[1])
    fold_size = length(permutation)÷n_folds

    #Get ranges for the cross validation folds. 
    folds = [((n-1)*fold_size + 1):(n*fold_size) for n in range(1, n_folds)]

    #Initialize class probabilisites
    PSMs[:,:prob] = zeros(Float64, size(PSMs)[1])
    #XGBoost model
    bst = ""
    for test_fold_idx in range(1, n_folds)
        train_fold_idxs = vcat([folds[fold] for fold in range(1, length(folds)) if fold != test_fold_idx]...)
        println(train_fold_idxs)
        train_fold_idxs = train_fold_idxs[rand(1:length(train_fold_idxs), min(max_train_size, length(train_fold_idxs)))]
        println(train_fold_idxs)
        return
        train_features = X[train_fold_idxs,:]
        train_classes = X_labels[train_fold_idxs,1]

        #Train a model on the n-1 training folds. Then apply it to get class probabilities for the test-fold. 
        bst = xgboost((train_features, train_classes), 
                        num_round=num_round, 
                        colsample_bytree = colsample_bytree, 
                        gamma = gamma, 
                        max_depth=max_depth, 
                        eta = eta, 
                        min_child_weight = min_child_weight, 
                        subsample = subsample, 
                        objective="binary:logistic",
                        #watchlist=(;)
                        )
        ŷ = XGBoost.predict(bst, X[folds[test_fold_idx],:])
        PSMs[folds[test_fold_idx],:prob] = (1 .- ŷ)
    end
    bst.feature_names = [string(x) for x in features]
    return bst
end

function rankPSMs!(PSMs::DataFrame, features::Vector{Symbol}; n_folds::Int = 3, colsample_bytree::Float64 = 0.5, num_round::Int = 25, eta::Float64 = 0.15, min_child_weight::Int = 1, subsample::Float64 = 0.5, gamma::Int = 0, max_depth::Int = 10, max_train_size::Int = 1e6, print_importance::Bool = false)
   
    #[:hyperscore,:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all, :spectral_contrast_matched,:RT_error,:scribe_score,:y_ladder,:b_ladder,:RT,:diff_hyper,:median_ions,:n_obs,:diff_scribe,:charge,:city_block,:matched_ratio,:weight,:intensity,:count,:SN]
    X = Matrix(PSMs[:,features])
    X_labels = PSMs[:, :decoy]

    #Train test split such that PSMs corresponding to the same peptide are always within the same 
    #cross-validation fold

    #Initialize class probabilisites
    PSMs[:,:prob] = zeros(Float64, size(PSMs)[1])

    folds = ones(Int64, size(PSMs)[1])
    sort!(PSMs, :sequence);
    fold = 1
    fold_ids = Tuple([_ for _ in range(1, n_folds)])
    for i in range(2, length(folds))
        if PSMs[i,:sequence] == PSMs[i-1,:sequence]
            folds[i] = fold
        else
            fold = rand(fold_ids)
            folds[i] = fold
        end
    end
    #XGBoost model
    println(unique(folds))
    bst = ""
    for test_fold_idx in range(1, n_folds)

        train_fold_idxs = findall(x->x!=test_fold_idx, folds)
        test_fold_idxs = findall(x->x==test_fold_idx, folds)

        train_features = X[train_fold_idxs,:]
        train_classes = X_labels[train_fold_idxs,1]

        #train_features = X[folds.!=test_fold_idx,:]
        #train_classes = X_labels[folds.!=test_fold_idx,1]
        #Train a model on the n-1 training folds. Then apply it to get class probabilities for the test-fold. 
        bst = xgboost((train_features, train_classes), 
                        num_round=num_round, 
                        colsample_bytree = colsample_bytree, 
                        gamma = gamma, 
                        max_depth=max_depth, 
                        eta = eta, 
                        min_child_weight = min_child_weight, 
                        subsample = subsample, 
                        objective="binary:logistic",
                        watchlist=(;)
                        )
        ŷ = XGBoost.predict(bst, X[test_fold_idxs,:])
        PSMs[test_fold_idxs,:prob] = (1 .- ŷ)
    end
    bst.feature_names = [string(x) for x in features]
    return bst, folds
end
