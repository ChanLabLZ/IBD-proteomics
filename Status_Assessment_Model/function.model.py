
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

def train_logistic_regression(X_train, y_train, X_test, y_test,
                               selected_features, selected_indices=None,
                               lasso_coef=None, top_n=10,
                               train_names=None, test_names=None,
                               model_prefix='lr'):
    """
    LogisticRegressionCV

    Return:
        lr_model, metrics, final_features
    """
    print("=" * 60)
    print("Step2：LogisticRegressionCV")
    print("=" * 60)

    n_selected = len(selected_features)
    n_train = X_train.shape[0]

    if n_selected > top_n and selected_indices is not None and lasso_coef is not None:
        coef_abs = np.abs(lasso_coef[selected_indices])
        top_indices_in_selected = np.argsort(coef_abs)[-top_n:][::-1]
        final_features = [selected_features[i] for i in top_indices_in_selected]

        print(f"LASSO selected {n_selected} > {top_n}, just using {top_n}.")
        print(f"Final features: {final_features}")

        X_train_final = X_train[:, top_indices_in_selected]
        X_test_final = X_test[:, top_indices_in_selected]
    else:
        final_features = selected_features
        X_train_final = X_train
        X_test_final = X_test
        print(f"LASSO selected {n_selected} (<={top_n}), using all.")

    lr_model = LogisticRegressionCV(
        Cs=10,
        cv=5,
        penalty='l2',
        solver='lbfgs',
        max_iter=5000,
        scoring='roc_auc',
        random_state=42)

    # avoid samples are too low
    if n_train <= 20:
        loo = LeaveOneOut()
        cv_scores = cross_val_score(lr_model, X_train_final, y_train, cv=loo, scoring='roc_auc')
        print(f"LOOCV AUC: {cv_scores.mean():.4f} (+/- {cv_scores.std()*2:.4f})")
    else:
        cv_folds = min(n_train, 5)
        cv_scores = cross_val_score(lr_model, X_train_final, y_train, cv=cv_folds, scoring='roc_auc')
        print(f"{cv_folds}-fold CV AUC: {cv_scores.mean():.4f} (+/- {cv_scores.std()*2:.4f})")

    # traing the final model
    lr_model.fit(X_train_final, y_train)

    # print the model parameter
    print(f"\nLogisticRegressionCV best C: {lr_model.C_[0]:.4f}")
    print(f"C search range: [{lr_model.Cs_[0]:.4f}, {lr_model.Cs_[-1]:.4f}]")

    # the performance of traing
    y_train_pred = lr_model.predict(X_train_final)
    y_train_proba = lr_model.predict_proba(X_train_final)[:, 1]

    train_auc = roc_auc_score(y_train, y_train_proba)
    train_cm = confusion_matrix(y_train, y_train_pred)

    # saving the resuts to review
    train_pred_df = pd.DataFrame({
        'sample': train_names if train_names is not None else range(len(y_train)),
        'true_label': y_train.astype(int),
        'pred_label': y_train_pred.astype(int),
        'pred_proba': np.round(y_train_proba, 6)})
    train_pred_df.to_csv(f'train_predictions_{model_prefix}.csv', index=False)

    train_cm_df = pd.DataFrame({
        '': ['Actual 0', 'Actual 1'],
        'Pred 0': [train_cm[0, 0], train_cm[1, 0]],
        'Pred 1': [train_cm[0, 1], train_cm[1, 1]]})
    train_cm_df.to_csv(f'train_confusion_matrix_{model_prefix}.csv', index=False)

    print("\n" + "-" * 40)
    print("The performance of model:")
    print("-" * 40)
    print(f"  AUC:         {train_auc:.4f}")
    print(f"  Confusion matrix:    TN={train_cm[0,0]} FP={train_cm[0,1]} | "
          f"FN={train_cm[1,0]} TP={train_cm[1,1]}")
    print(f"Saving: train_predictions_{model_prefix}.csv, train_confusion_matrix_{model_prefix}.csv")

    # the performance of test
    y_test_pred = lr_model.predict(X_test_final)
    y_test_proba = lr_model.predict_proba(X_test_final)[:, 1]

    metrics = {
        'train_auc': train_auc,
        'accuracy':  accuracy_score(y_test, y_test_pred),
        'precision': precision_score(y_test, y_test_pred, zero_division=0),
        'recall':    recall_score(y_test, y_test_pred, zero_division=0),
        'f1':        f1_score(y_test, y_test_pred, zero_division=0),
        'auc':       roc_auc_score(y_test, y_test_proba)}

    test_cm = confusion_matrix(y_test, y_test_pred)

    # saving the test resuts to review
    test_pred_df = pd.DataFrame({
        'sample': test_names if test_names is not None else range(len(y_test)),
        'true_label': y_test_pred.astype(int),
        'pred_label': y_test_proba.astype(int),
        'pred_proba': np.round(y_test_proba, 6)})
    test_pred_df.to_csv(f'test_predictions_{model_prefix}.csv', index=False)

    test_cm_df = pd.DataFrame({
        '': ['Actual 0', 'Actual 1'],
        'Pred 0': [test_cm[0, 0], test_cm[1, 0]],
        'Pred 1': [test_cm[0, 1], test_cm[1, 1]]})
    test_cm_df.to_csv(f'test_confusion_matrix_{model_prefix}.csv', index=False)

    print("\n" + "-" * 40)
    print("Performance of test set:")
    print("-" * 40)
    for k, v in metrics.items():
        print(f"  {k:>12s}: {v:.4f}")

    print(f"  Confusion matrix:    TN={test_cm[0,0]} FP={test_cm[0,1]} | "
          f"FN={test_cm[1,0]} TP={test_cm[1,1]}")
    print(f"Saving: test_predictions_{model_prefix}.csv, test_confusion_matrix_{model_prefix}.csv")

    # Coef
    coef_df = pd.DataFrame({
        'feature': final_features,
        'coefficient': lr_model.coef_[0]
    }).sort_values('coefficient', key=abs, ascending=False)
    print(f"\nLR coef (top 10):")
    print(coef_df.head(10).to_string(index=False))
    coef_df.to_csv(f'lr_coefficients_{model_prefix}.csv', index=False)

    print()
    return (lr_model, metrics, final_features,
            y_train_proba, y_test_proba,
            y_train_pred, y_test_pred,
            X_train_final, X_test_final)


def load_model(model_path):
    """loading bundle"""
    print(f"Loading model from: {model_path}")
    bundle = joblib.load(model_path)
    model_key = 'rf_model' if 'rf_model' in bundle else 'lr_model'
    model = bundle[model_key]
    scaler = bundle.get('scaler')
    model_name = 'RF' if model_key == 'rf_model' else 'LR'
    print(f"Model type: {model_name}")
    print(f"Features used: {len(bundle['final_features'])}")
    print(f"Feature list: {bundle['final_features']}")
    return bundle, model, scaler, model_name


def predict_and_evaluate(model, X_val, y_true, sample_ids, model_name):

    print(f"\n{'='*60}")
    print(f"{model_name} Validation Results")
    print(f"{'='*60}")

    y_pred = model.predict(X_val)
    y_proba = model.predict_proba(X_val)[:, 1]

    # metrics
    metrics = {
        'accuracy':  accuracy_score(y_true, y_pred),
        'precision': precision_score(y_true, y_pred, zero_division=0),
        'recall':    recall_score(y_true, y_pred, zero_division=0),
        'f1':        f1_score(y_true, y_pred, zero_division=0),
    }
    try:
        metrics['auc'] = roc_auc_score(y_true, y_proba)
    except ValueError:
        metrics['auc'] = float('nan')

    cm = confusion_matrix(y_true, y_pred)

    for k, v in metrics.items():
        print(f"  {k:>12s}: {v:.4f}")
    print(f"  Confusion matrix:    TN={cm[0,0]} FP={cm[0,1]} | FN={cm[1,0]} TP={cm[1,1]}")

    # saving validation results
    pred_df = pd.DataFrame({
        'sample': sample_ids,
        'true_label': y_true.astype(int),
        'pred_label': y_pred.astype(int),
        'pred_proba': np.round(y_proba, 6),
    })
    pred_df.to_csv('validation_predictions.csv', index=False)
    print("Saving: validation_predictions.csv")

    # Saving confusion
    cm_df = pd.DataFrame({
        '': ['Actual 0', 'Actual 1'],
        'Pred 0': [cm[0, 0], cm[1, 0]],
        'Pred 1': [cm[0, 1], cm[1, 1]],
    })
    cm_df.to_csv('validation_confusion_matrix.csv', index=False)
    print("Saving: validation_confusion_matrix.csv")

    from function.plot import plot_val_roc, plot_val_cm
    
    # ROC plot
    plot_val_roc(y_true, y_proba, model_name)

    # plot cpmfusion
    plot_val_cm(y_true, y_pred, model_name)

    return metrics
