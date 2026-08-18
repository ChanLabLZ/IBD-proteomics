
import numpy as np
import pandas as pd
import joblib
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, confusion_matrix
)
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
import warnings
warnings.filterwarnings('ignore')


from function.model import load_model, predict_and_evaluate
from function.plot import plot_val_roc, plot_val_cm



MODEL_PATH = '../final_model_lr_aTNF.pkl'
MATRIX_PATH = '../..//protein_Samplematrix_imputeNA1.xlsx'
SAMPLE_INFO_PATH = '../../XHOM_PM_sampleInfor.xlsx'

DRUG = 'anti-TNF'
LABEL_COL = 'Group'  # R / NR
SAMPLE_ID_COL = 'Type'


def main():
    print("\n" + "=" * 60)
    print("Model Validation")
    print("=" * 60 + "\n")

    # ---- 1. Loading model ----
    bundle, model, scaler, model_name = load_model(MODEL_PATH)

    # ---- 2. Loading data ----
    print("Loading new data...")

    # matrix
    discovery_set = pd.read_excel(MATRIX_PATH, index_col='protein')
    discovery_set = discovery_set.T

    DEPs=pd.read_csv('./a1_PostvsPre_RNR_deg2.txt', sep='\t', header=0)
    DEPs=DEPs[DEPs['type']!='NOT']
    discovery_set=discovery_set.loc[:,DEPs['protein']]

    # sample information
    sampleInfor = pd.read_excel(SAMPLE_INFO_PATH)
    sampleInfor = sampleInfor.loc[
        sampleInfor['Group'].isin(['R', 'NR']) &
        (sampleInfor['Description'] == DRUG)]
    print(sampleInfor['Group'].value_counts(), "\n")

    sampleInfor['state'] = sampleInfor.apply(lambda r: 0 if r['Group'] == 'NR' else 1, axis=1)
    print(sampleInfor['state'].value_counts(), "\n")

    discovery_set = discovery_set.loc[sampleInfor[SAMPLE_ID_COL], :]
    feature_names = list(discovery_set.columns)

    X_raw = discovery_set
    y_val = np.array(sampleInfor['state'], dtype=np.float64)
    sample_ids = sampleInfor[SAMPLE_ID_COL].values

    print(f"Sample counts: {len(y_val)}")

    # ---- 3. Align features and normalize ----

    if scaler is not None:
        X_val = pd.DataFrame(
            scaler.transform(X_raw),
            index=X_raw.index,
            columns=X_raw.columns)
    else:
        X_val = X_raw
    
    final_features = bundle['final_features']
    missing = [f for f in final_features if f not in X_raw.columns]
    if missing:
        print(f"Warnings: these features: {missing} are missing in the matrix.")

    available = [f for f in final_features if f in X_raw.columns]
    if len(available) != len(final_features):
        print(f"Warnings: need {len(final_features)} features，only {len(available)} could be used.")

    X_sel = X_val[available]

    print(f"Final features: {len(available)}")

    # ---- 4. Prediction and validation ----
    metrics = predict_and_evaluate(model, X_sel, y_val, sample_ids, model_name)

    print(f"\n{'='*60}")
    print("Validation finishing!")
    print("="*60)


if __name__ == '__main__':
    main()
