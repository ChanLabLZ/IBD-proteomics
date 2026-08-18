
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.model_selection import train_test_split, cross_val_score, LeaveOneOut
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LassoCV, LogisticRegressionCV, ElasticNet
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RandomizedSearchCV
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, confusion_matrix
)
import matplotlib.pyplot as plt
import joblib
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

from function.selection import lasso_feature_selection
from function.model import train_logistic_regression
from function.plot import plot_feature_importance,plot_lr_coefficient_bar,plot_roc_curve,plot_confusion_matrices

def load_data(matrix_path, sample_path, feature_path=None, drug='anti-TNF', sel_group='protein'):
    if feature_path:
        DEPs = pd.read_csv(feature_path, sep='\t', header=0)
        DEPs = DEPs[DEPs['type'] != 'NOT']
    else:
        print("No feature files！")

    discovery_set = pd.read_excel(matrix_path, index_col='protein')
    discovery_set = discovery_set.T

    if feature_path and sel_group in DEPs.columns:
        discovery_set = discovery_set.loc[:, DEPs[sel_group]]

    sampleInfor = pd.read_excel(sample_path)
    sampleInfor = sampleInfor.loc[
        sampleInfor['Group'].isin(['NR', 'R']) &
        (sampleInfor['Description'] == drug),
    ]
    print(sampleInfor['Group'].value_counts())
    sampleInfor['state'] = sampleInfor.apply(
        lambda r: 0 if r['Group'] == 'NR' else 1, axis=1
    )
    print(sampleInfor['state'].value_counts())
    return discovery_set, sampleInfor


def shap_analysis(lr_model, X_train, X_test, feature_names, prefix='lr'):
    """
    SHAP analysis (LR)：LinearExplainer explain Logistic Regression prediction。
    """
    import shap

    print("=" * 60)
    print("SHAP analysis of LR")
    print("=" * 60)

    X_train_df = pd.DataFrame(X_train, columns=feature_names)
    X_test_df = pd.DataFrame(X_test, columns=feature_names)

    explainer = shap.LinearExplainer(lr_model, X_train_df)
    shap_train = explainer(X_train_df)
    shap_test = explainer(X_test_df)

    # ---------- training ----------
    fig, ax = plt.subplots(figsize=(8, max(4, len(feature_names) * 0.3)))
    shap.plots.bar(shap_train, max_display=len(feature_names), show=False)
    ax.set_title(f'SHAP Feature Importance (Train, n={X_train.shape[0]})')
    plt.tight_layout()
    plt.savefig(f'shap_bar_train_{prefix}.pdf', bbox_inches='tight')
    print(f"Saving: shap_bar_train_{prefix}.pdf")
    plt.close()

    fig, ax = plt.subplots(figsize=(8, max(4, len(feature_names) * 0.3)))
    shap.plots.beeswarm(shap_train, max_display=len(feature_names), show=False)
    ax.set_title(f'SHAP Beeswarm (Train, n={X_train.shape[0]})')
    plt.tight_layout()
    plt.savefig(f'shap_beeswarm_train_{prefix}.pdf', bbox_inches='tight')
    print(f"Saving: shap_beeswarm_train_{prefix}.pdf")
    plt.close()

    shap_train_df = pd.DataFrame(
        shap_train.values,
        columns=[f'{f}_shap' for f in feature_names]
    )
    shap_train_df.insert(0, 'expected_value', shap_train.base_values)
    shap_train_df.to_csv(f'shap_values_train_{prefix}.csv', index=False)
    print(f"Saving: shap_values_train_{prefix}.csv")

    # ---------- Test ----------
    fig, ax = plt.subplots(figsize=(8, max(4, len(feature_names) * 0.3)))
    shap.plots.bar(shap_test, max_display=len(feature_names), show=False)
    ax.set_title(f'SHAP Feature Importance (Test, n={X_test.shape[0]})')
    plt.tight_layout()
    plt.savefig(f'shap_bar_test_{prefix}.pdf', bbox_inches='tight')
    print(f"Saving: shap_bar_test_{prefix}.pdf")
    plt.close()

    fig, ax = plt.subplots(figsize=(8, max(4, len(feature_names) * 0.3)))
    shap.plots.beeswarm(shap_test, max_display=len(feature_names), show=False)
    ax.set_title(f'SHAP Beeswarm (Test, n={X_test.shape[0]})')
    plt.tight_layout()
    plt.savefig(f'shap_beeswarm_test_{prefix}.pdf', bbox_inches='tight')
    print(f"Saving: shap_beeswarm_test_{prefix}.pdf")
    plt.close()

    shap_test_df = pd.DataFrame(
        shap_test.values,
        columns=[f'{f}_shap' for f in feature_names]
    )
    shap_test_df.insert(0, 'expected_value', shap_test.base_values)
    shap_test_df.to_csv(f'shap_values_test_{prefix}.csv', index=False)
    print(f"Saving: shap_values_test_{prefix}.csv")

    print()
    return shap_train, shap_test


def main():
    print("LASSO + LogisticRegressionCV")

    # ---- 1. loading data ----
    discovery_set = pd.read_excel('../../protein_Samplematrix_imputeNA1.xlsx',index_col='protein')
    discovery_set = discovery_set.T

    DEPs=pd.read_csv('a1_PostvsPre_RNR_deg2.txt', sep='\t', header=0)
    DEPs=DEPs[DEPs['type']!='NOT']

    discovery_set=discovery_set.loc[:,DEPs['protein']]

    sampleInfor = pd.read_excel('../../XHOM_PM_sampleInfor.xlsx')
    sampleInfor = sampleInfor.loc[
        sampleInfor['Group'].isin(['R','NR']) &
        (sampleInfor['Description'] == 'anti-TNF') &
        (sampleInfor['Group_prepost'] != 'post')]
    print(sampleInfor['Group'].value_counts(), "\n")

    sampleInfor['state'] = sampleInfor.apply(lambda r: 0 if r['Group'] == 'NR' else 1, axis=1)
    print(sampleInfor['state'].value_counts(), "\n")

    # ---- 2. training/test set ----
    X_train_idx, X_val_idx, y_train_idx, y_val_idx = train_test_split(sampleInfor, sampleInfor['state'],
        test_size=0.3, random_state=260520, stratify=sampleInfor['state'])

    discovery_set = discovery_set.loc[sampleInfor['Type'], :]
    feature_names = list(discovery_set.columns)

    X_train = discovery_set.loc[X_train_idx['Type'], :]
    y_train = np.array(y_train_idx, dtype=np.float64)
    X_test = discovery_set.loc[X_val_idx['Type'], :]
    y_test = np.array(y_val_idx, dtype=np.float64)

    print(f"All features: {len(feature_names)}")
    print(f"Trainging: {X_train.shape[0]} 样本")
    print(f"Test: {X_test.shape[0]} 样本")
    print(f"train NR(0)/R(1) = {sum(y_train==0)}/{sum(y_train==1)}  "
          f"test NR(0)/R(1) = {sum(y_test==0)}/{sum(y_test==1)}\n")

    # ---- 3. pre-filter ----
    ## not applying, we have used the DEPs
    filtered_features = feature_names
    X_train_filtered = X_train[filtered_features]
    X_test_filtered = X_test[filtered_features]

    # ---- 4. feature selection by LASSO ----
    (X_train_selected, X_test_selected, selected_features, lasso_model, scaler,
    selected_indices, lasso_coef) = lasso_feature_selection(X_train_filtered, y_train, X_test_filtered, filtered_features)

    if len(selected_features) == 0:
        print("Warnings: 0 features were selected by LASSO! Will using all features.")
        selected_features = filtered_features
        X_train_selected = scaler.transform(X_train_filtered)
        X_test_selected = scaler.transform(X_test_filtered)
        selected_indices = np.arange(len(filtered_features))
        lasso_coef = np.ones(len(filtered_features))

    # ---- 5. Train LR and validation ----
    train_names = X_train.index.values
    test_names = X_test.index.values

    (lr_model, lr_metrics, lr_final_features, lr_y_train_proba, lr_y_test_proba,
    lr_y_train_pred, lr_y_test_pred, lr_X_train_final, lr_X_test_final) = train_logistic_regression(
        X_train_selected, y_train,
        X_test_selected, y_test,
        selected_features,
        selected_indices=selected_indices,
        lasso_coef=lasso_coef,
        top_n=10,
        train_names=train_names,
        test_names=test_names,
        model_prefix='lr')

    # ---- 6. Visualization ----
    print("Ploting...")

    plot_feature_importance(lasso_model, lr_model, filtered_features, lr_final_features, prefix='lr')
    plot_lr_coefficient_bar(lr_model, lr_final_features, prefix='lr')
    shap_analysis(lr_model, lr_X_train_final, lr_X_test_final, lr_final_features, prefix='lr')
    plot_roc_curve(y_train, lr_y_train_proba, y_test, lr_y_test_proba, prefix='lr')
    plot_confusion_matrices(y_train, lr_y_train_pred, y_test, lr_y_test_pred, prefix='lr')


    # ---- 7. Summary ----
    print("\n" + "=" * 60)
    print("Process summary")
    print(f"All features:     {len(feature_names)}")
    print(f"Pre-filter:       {len(filtered_features)}")
    print(f"Slected by LASSO:     {len(selected_features)}")
    print(f"Using in LR:   {len(lr_final_features)}")

    print("The performance of LR")
    print(f"  {'Index':<12s} {'LR':>8s}")
    print(f"  {'-'*30}")
    for k in ['train_auc', 'auc', 'accuracy', 'precision', 'recall', 'f1']:
        print(f"  {k:<12s} {lr_metrics[k]:>8.4f}")
    print(f"{'='*40}")
    print("=" * 60)

    # saving model
    print("\n"+"Saving model...")

    lr_bundle = {
        'lr_model': lr_model,
        'scaler': scaler,
        'final_features': lr_final_features,
        'selected_features': selected_features,
        'filtered_features': filtered_features,
        'feature_names': feature_names,
        'metrics': lr_metrics}
    joblib.dump(lr_bundle, 'final_model_lr_aTNF.pkl')
    print("Saving: final_model_lr_aTNF.pkl")

    return lasso_model, lr_model, lr_metrics, lr_final_features

if __name__ == "__main__":
    lasso_model, lr_model, lr_metrics, lr_final_features = main()
