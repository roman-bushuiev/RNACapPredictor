import numpy as np
import pandas as pd
from sklearn.neighbors import KNeighborsClassifier


# Function to create feature vector for each sample
def create_features(df, all_rt_names, include_insdel):
    features = []
    labels = []
    caps = []
    experiments = []
    if include_insdel:
        nucls = ['A%_INSDEL', 'C%_INSDEL', 'G%_INSDEL', 'T%_INSDEL', 'INS%_INSDEL', 'DEL%_INSDEL']
    else:
        nucls = ['A%', 'C%', 'G%', 'T%']
    
    for nucl in nucls:
        if nucl not in df.columns:
            raise ValueError(f"Column {nucl} not found in DataFrame")
    
    # Group by cap type and experiment
    for (cap_type, experiment), group in df.groupby(['cap', 'experiment']):
        # Initialize feature vector with zeros
        feature_vec = np.zeros(len(all_rt_names) * len(nucls))

        # For each RT in this group, add its ACGT percentages
        for _, row in group.iterrows():
            rt_idx = np.where(all_rt_names == row['RT'])[0][0]
            base_idx = rt_idx * len(nucls)
            feature_vec[base_idx:base_idx+len(nucls)] = row[nucls]
        
        features.append(feature_vec)
        labels.append(f"{cap_type} ({experiment})")
        caps.append(cap_type)
        experiments.append(experiment)

    return np.array(features), np.array(labels), np.array(caps), np.array(experiments)


# Predict using k-NN with cosine similarity and masked training set depending on the RTs present in the test sample
def predict(X_test_sample, X_train, y_train):
    knn = KNeighborsClassifier(n_neighbors=X_train.shape[0], metric='cosine')
    mask = X_test_sample != 0
    X_train_masked = X_train.copy()
    X_train_masked[:, ~mask] = 0  # Use only training RTs that are present in the test sample
    knn.fit(X_train_masked[:, mask], y_train)
    
    # Get distances and indices of 5 nearest neighbors
    distances, indices = knn.kneighbors(X_test_sample[mask].reshape(1, -1))
    
    # Convert distances to similarities (1 - distance)
    similarities = 1 - distances[0]
    
    # Get the corresponding labels
    neighbor_labels = y_train[indices[0]]
    
    return list(zip(neighbor_labels, similarities))


def predict_cap(df_train, df_test, show_true_cap=False, include_insdel=False, print_top_k=50):
    # Get unique RT names from all datasets
    all_rt_names = df_train['RT'].unique()

    # Create features and labels for each dataset
    X_train, y_train, caps_train, experiments_train = create_features(df_train, all_rt_names, include_insdel)
    X_test, y_test, caps_test, experiments_test = create_features(df_test, all_rt_names, include_insdel)

    # Make predictions
    test_predictions = [predict(x, X_train, y_train) for x in X_test]

    # Create a list to store results for DataFrame
    results = []
    
    # Print predictions with similarities and build DataFrame
    for i, (true, preds) in enumerate(zip(y_test, test_predictions)):
        used_rts = df_test[df_test['experiment'] == experiments_test[i]]['RT']
        mean_reads = df_test[df_test['experiment'] == experiments_test[i]]['num_reads_ACGT'].mean()
        
        print(f"Experiment: {experiments_test[i]}")
        if show_true_cap:
            print(f"True cap: {caps_test[i]}")
        print(f"{len(used_rts)} RTs considered for prediction({used_rts.tolist()}) with mean "
              f"number of reads {mean_reads}")
        
        # Store results for each prediction
        result_dict = {
            'experiment': experiments_test[i],
            'true_cap': caps_test[i] if show_true_cap else None,
            'num_rts': len(used_rts),
            'used_rts': used_rts.tolist(),
            'mean_reads': mean_reads
        }
        
        for k, (pred, sim) in enumerate(preds[:print_top_k]):
            print(f"Top-{k+1} prediction: {pred:8} with similarity {sim:.3f}")
            result_dict[f'prediction_{k+1}'] = pred
            result_dict[f'similarity_{k+1}'] = sim
            
        results.append(result_dict)
        print("\n")
    
    # Create DataFrame from results
    results_df = pd.DataFrame(results)
    return results_df


def main():
    raise NotImplementedError("Not implemented yet.")
    # df_train = pd.read_csv('data/train.csv')
    # df_test = pd.read_csv('data/test.csv')
    # predict_cap(df_train, df_test)


if __name__ == "__main__":
    main()
