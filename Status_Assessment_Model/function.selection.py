
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


def lasso_feature_selection(X_train, y_train, X_test, feature_names):
    """
    LASSO selection

    Return:
        X_train_selected: (n_train, n_selected) 
        X_test_selected:  (n_test,  n_selected)
        selected_features: feature name
        lasso_cv: LASSO model
        scaler: StandardScaler
        selected_indices: feature index
        coef: LASSO coef
    """

    print("=" * 60)
    print("Step1：feature selection by LASSO")
    print("=" * 60)

    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    n_train = X_train_scaled.shape[0]
    cv_folds = min(n_train, 5)  # if the samples is too low, using LOOCV，else 5-fold
    if n_train <= 10:
        cv_strategy = LeaveOneOut()
    else:
        cv_strategy = cv_folds

    print(f"Samples: {n_train}, CV strategy: {'LOOCV' if n_train <= 10 else f'{cv_folds}-fold'}")

    lasso_cv = LassoCV(
        cv=cv_strategy,
        random_state=42,
        max_iter=20000,
        n_alphas=200,
        eps=1e-4,
        tol=1e-4)

    lasso_cv.fit(X_train_scaled, y_train)

    print(f"Best alpha value: {lasso_cv.alpha_:.6f}")

    coef = lasso_cv.coef_
    non_zero_mask = coef != 0
    selected_indices = np.where(non_zero_mask)[0]
    selected_features = [feature_names[i] for i in selected_indices]

    print(f"LASSO selected features: {len(selected_features)}/{len(feature_names)}")
    if len(selected_features) > 0:
        coef_abs = np.abs(coef[selected_indices])
        top5 = np.argsort(coef_abs)[-5:][::-1]
        print("The top 5 features by LASSO:")
        for idx in top5:
            fi = selected_indices[idx]
            print(f"  {feature_names[fi]}: coeff={coef[fi]:.4f}")

    X_train_selected = X_train_scaled[:, non_zero_mask]
    X_test_selected = X_test_scaled[:, non_zero_mask]

    print()
    return X_train_selected, X_test_selected, selected_features, lasso_cv, scaler, selected_indices, coef


