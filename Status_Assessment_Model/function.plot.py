
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



def plot_feature_importance(lasso_model, lr_model, feature_names, selected_features, prefix='lr'):
    """The importance of features and CSV"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    lasso_coef = lasso_model.coef_
    non_zero_idx = np.where(lasso_coef != 0)[0]

    if len(non_zero_idx) > 0:
        sorted_idx = non_zero_idx[np.argsort(np.abs(lasso_coef[non_zero_idx]))]
        axes[0].barh(range(len(sorted_idx)), lasso_coef[sorted_idx])
        axes[0].set_yticks(range(len(sorted_idx)))
        axes[0].set_yticklabels([feature_names[i] for i in sorted_idx], fontsize=6)
        axes[0].set_xlabel('Coefficient')
        axes[0].set_title(f'LASSO ({len(non_zero_idx)} features)')
        axes[0].axvline(x=0, color='black', linestyle='--', linewidth=0.5)

        pd.DataFrame({
            'feature': [feature_names[i] for i in sorted_idx],
            'coefficient': lasso_coef[sorted_idx],
            'coefficient_abs': np.abs(lasso_coef[sorted_idx])
        }).to_csv(f'lasso_coefficients_{prefix}.csv', index=False)
        print(f"Saving: lasso_coefficients_{prefix}.csv")

    lr_coef = lr_model.coef_[0]
    sorted_idx = np.argsort(np.abs(lr_coef))
    axes[1].barh(range(len(lr_coef)), lr_coef[sorted_idx])
    axes[1].set_yticks(range(len(lr_coef)))
    axes[1].set_yticklabels([selected_features[i] for i in sorted_idx], fontsize=7)
    axes[1].set_xlabel('Coefficient')
    axes[1].set_title('Logistic Regression')
    axes[1].axvline(x=0, color='black', linestyle='--', linewidth=0.5)

    pd.DataFrame({
        'feature': [selected_features[i] for i in sorted_idx],
        'coefficient': lr_coef[sorted_idx],
        'coefficient_abs': np.abs(lr_coef[sorted_idx])
    }).to_csv(f'lr_coefficients_{prefix}.csv', index=False)
    print(f"Saving: lr_coefficients_{prefix}.csv")

    plt.tight_layout()
    plt.savefig(f'feature_importance_{prefix}.pdf', bbox_inches='tight')
    print(f"Saving: feature_importance_{prefix}.pdf")
    plt.close()


def plot_lr_coefficient_bar(lr_model, feature_names, prefix='lr'):
    """obs value of LR coef bar plot and CSV"""
    coef = lr_model.coef_[0]
    coef_abs = np.abs(coef)
    sorted_idx = np.argsort(coef_abs)[::-1]  # ording

    n = len(coef)
    fig, ax = plt.subplots(figsize=(8, max(4, n * 0.35)))

    colors = ['#d62728' if coef[sorted_idx][i] < 0 else '#2ca02c'
              for i in range(n)]
    ax.barh(range(n), coef_abs[sorted_idx], color=colors)
    ax.set_yticks(range(n))
    ax.set_yticklabels([feature_names[i] for i in sorted_idx])
    ax.invert_yaxis()
    ax.set_xlabel('|Coefficient|')
    ax.set_title(f'Logistic Regression Feature Importance ({n} features)')

    # figure legend
    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(color='#2ca02c', label='Positive'),
        Patch(color='#d62728', label='Negative'),
    ], loc='lower right')

    plt.tight_layout()
    plt.savefig(f'lr_coefficient_bar_{prefix}.pdf', bbox_inches='tight')
    print(f"Saving: lr_coefficient_bar_{prefix}.pdf")

    # Saving CSV
    pd.DataFrame({
        'feature': [feature_names[i] for i in sorted_idx],
        'coefficient': coef[sorted_idx],
        'abs_coefficient': coef_abs[sorted_idx],
    }).to_csv(f'lr_coefficient_bar_{prefix}.csv', index=False)
    print(f"Saving: lr_coefficient_bar_{prefix}.csv")
    plt.close()


def plot_roc_curve(y_train, y_train_proba, y_test, y_test_proba, prefix='lr'):
    """ROC of training and test sets and CSV"""
    from sklearn.metrics import roc_curve as roc_curve_fn

    # ---- training set ----
    fpr_train, tpr_train, thr_train = roc_curve_fn(y_train, y_train_proba)
    auc_train = roc_auc_score(y_train, y_train_proba)

    pd.DataFrame({
        'fpr': fpr_train, 'tpr': tpr_train,
        'threshold': np.round(thr_train, 6), 'auc': auc_train, 'set': 'train'
    }).to_csv(f'roc_curve_train_{prefix}.csv', index=False)
    print(f"Saving: roc_curve_train_{prefix}.csv")

    # ---- test ----
    fpr_test, tpr_test, thr_test = roc_curve_fn(y_test, y_test_proba)
    auc_test = roc_auc_score(y_test, y_test_proba)

    pd.DataFrame({
        'fpr': fpr_test, 'tpr': tpr_test,
        'threshold': np.round(thr_test, 6), 'auc': auc_test, 'set': 'test'
    }).to_csv(f'roc_curve_test_{prefix}.csv', index=False)
    print(f"Saving: roc_curve_test_{prefix}.csv")

    # ---- merge ----
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.plot(fpr_train, tpr_train, color='steelblue', lw=2,
            label=f'Train (AUC = {auc_train:.4f})')
    ax.plot(fpr_test, tpr_test, color='darkorange', lw=2,
            label=f'Test  (AUC = {auc_test:.4f})')
    ax.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    ax.set_xlim([0.0, 1.0]); ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title(f'ROC Curve  |  Train AUC={auc_train:.4f}  Test AUC={auc_test:.4f}')
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'roc_curve_{prefix}.pdf', bbox_inches='tight')
    print(f"Saving: roc_curve_{prefix}.pdf")
    plt.close()

def plot_confusion_matrices(y_train, y_train_pred, y_test, y_test_pred, prefix='lr'):
    from sklearn.metrics import confusion_matrix as cm_fn

    cm_train = cm_fn(y_train, y_train_pred)
    cm_test = cm_fn(y_test, y_test_pred)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax, cm, title, set_name in [
        (axes[0], cm_train, 'Train', 'train'),
        (axes[1], cm_test,  'Test',  'test'),
    ]:
        im = ax.imshow(cm, cmap='Blues', interpolation='nearest')
        ax.set_title(f'{title} Confusion Matrix')
        ax.set_xlabel('Predicted'); ax.set_ylabel('Actual')
        ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
        ax.set_xticklabels(['0', '1']); ax.set_yticklabels(['0', '1'])

        for i in range(2):
            for j in range(2):
                ax.text(j, i, str(cm[i, j]),
                        ha='center', va='center',
                        fontsize=20, fontweight='bold',
                        color='white' if cm[i, j] > cm.max() / 2 else 'black')

        # Saving CSV
        pd.DataFrame({
            '': ['Actual 0', 'Actual 1'],
            'Pred 0': [cm[0, 0], cm[1, 0]],
            'Pred 1': [cm[0, 1], cm[1, 1]],
        }).to_csv(f'confusion_matrix_{set_name}_{prefix}.csv', index=False)
        print(f"Saving: confusion_matrix_{set_name}_{prefix}.csv")

    plt.tight_layout()
    plt.savefig(f'confusion_matrices_{prefix}.pdf', bbox_inches='tight')
    print(f"Saving: confusion_matrices_{prefix}.pdf")
    plt.close()


def plot_val_roc(y_val, y_val_proba, model_name):
    fpr, tpr, _ = roc_curve(y_val, y_val_proba)
    auc = roc_auc_score(y_val, y_val_proba)

    pd.DataFrame({
        'fpr': fpr, 'tpr': tpr,
        'threshold': np.round(_, 6), 'auc': auc, 'set': 'validation'
    }).to_csv('roc_curve_val.csv', index=False)

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.plot(fpr, tpr, color='darkorange', lw=2, label=f'{model_name} Val (AUC = {auc:.4f})')
    ax.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    ax.set_xlim([0.0, 1.0]); ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title(f'{model_name} Validation ROC Curve')
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('roc_curve_val.pdf', bbox_inches='tight')
    plt.close()


def plot_val_cm(y_true, y_pred, model_name):
    cm = confusion_matrix(y_true, y_pred)
    fig, ax = plt.subplots(figsize=(5, 5))
    im = ax.imshow(cm, cmap='Blues', interpolation='nearest')
    ax.set_title(f'{model_name} Validation Confusion Matrix')
    ax.set_xlabel('Predicted'); ax.set_ylabel('Actual')
    ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
    ax.set_xticklabels(['0', '1']); ax.set_yticklabels(['0', '1'])
    for i in range(2):
        for j in range(2):
            ax.text(j, i, str(cm[i, j]),
                    ha='center', va='center',
                    fontsize=20, fontweight='bold',
                    color='white' if cm[i, j] > cm.max() / 2 else 'black')
    plt.tight_layout()
    plt.savefig('confusion_matrix_val.pdf', bbox_inches='tight')
    plt.close()

