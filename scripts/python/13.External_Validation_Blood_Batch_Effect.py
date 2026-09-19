
# EXTERNAL VALIDATION ON A SECOND, INDEPENDENT BLOOD COHORT (GSE10715)


import sys
import os
import pandas as pd
import numpy as np
import gzip
import joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, roc_auc_score, accuracy_score, precision_score, recall_score, f1_score, confusion_matrix
import matplotlib.pyplot as plt


if len(sys.argv) < 6:
    print("Error: Missing file arguments!")
    sys.exit(1)

importance_csv_file = sys.argv[1]
matrix_164191_file   = sys.argv[2]
gpl_platform_file    = sys.argv[3]
matrix_10715_file     = sys.argv[4]
output_dir            = sys.argv[5]




def bootstrap_ci_binary(y_true, y_prob, n_boot=2000, seed=123):
    rng = np.random.RandomState(seed)
    y_true = np.array(y_true)
    y_prob = np.array(y_prob)
    n = len(y_true)
    aucs, accs, sens, precs, f1s = [], [], [], [], []
    for _ in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        yt, yp = y_true[idx], y_prob[idx]
        if len(np.unique(yt)) < 2:
            continue
        pred = (yp >= 0.5).astype(int)
        aucs.append(roc_auc_score(yt, yp))
        accs.append(accuracy_score(yt, pred))
        sens.append(recall_score(yt, pred, zero_division=0))
        precs.append(precision_score(yt, pred, zero_division=0))
        f1s.append(f1_score(yt, pred, zero_division=0))
    def ci(x): return np.percentile(x, [2.5, 97.5])
    return {"AUC_CI": ci(aucs), "Accuracy_CI": ci(accs),
            "Sensitivity_CI": ci(sens), "Precision_CI": ci(precs), "F1_CI": ci(f1s)}


def load_data():
    print("Loading GSE164191 (Blood Data)...")
    file_path = matrix_164191_file
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Dataset not found at {file_path}")

    df = pd.read_csv(file_path, sep='\t', comment='!', index_col=0)

    with gzip.open(file_path, 'rt') as f:
        for line in f:
            if line.startswith('!Sample_title'):
                titles = line.strip().split('\t')[1:]
                break

    y = np.array([1 if any(x in t.lower() for x in ['crc', 'cancer', 'colorectal']) else 0 for t in titles])
    print(f"Samples: {len(y)}, Cancer: {sum(y)}, Healthy: {len(y)-sum(y)}")
    return df, y


def get_gene_mapping():
    file_path = gpl_platform_file
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"GPL570 mapping file not found at {file_path}")

    gpl = pd.read_csv(file_path, sep='\t', comment='#', low_memory=False)
    mapping = {str(row['ID']): str(row['Gene Symbol']).split(' /// ') for _, row in gpl.iterrows()}
    return mapping


df_importance = pd.read_csv(importance_csv_file)
target_genes = df_importance["Gene"].tolist()


def select_features(X_train, y_train, mapping):
    selected = []
    cancer_idx = np.where(y_train == 1)[0]
    healthy_idx = np.where(y_train == 0)[0]

    for probe in X_train.index:
        symbols = mapping.get(probe, [])
        if not any(s in target_genes for s in symbols):
            continue

        c = X_train.loc[probe].iloc[cancer_idx].values.astype(float)
        h = X_train.loc[probe].iloc[healthy_idx].values.astype(float)

        if len(c) < 2 or len(h) < 2:
            continue

        fc = np.mean(c) / (np.mean(h) + 1e-6)
        from scipy import stats
        p = stats.ttest_ind(c, h, equal_var=False).pvalue

        if (fc > 1.1 or fc < 0.9) and p < 0.05:
            selected.append(probe)
    return selected



# Train ONE final, locked blood model on the FULL GSE164191 dataset
# GSE10715 is a fully independent cohort 


df, y = load_data() 
mapping = get_gene_mapping() 

final_selected_probes = select_features(df, y, mapping) 
print(f"Final locked blood panel: {len(final_selected_probes)} probes")

final_model = RandomForestClassifier(n_estimators=100, random_state=42)
final_model.fit(df.loc[final_selected_probes].T, y)

# Save this as the frozen
joblib.dump(final_model, os.path.join(output_dir, "Final_Locked_Blood_Model.pkl"))
pd.Series(final_selected_probes).to_csv(os.path.join(output_dir, "Final_Locked_Blood_Probes.csv"), index=False)





def load_gse10715():
    file_path = matrix_10715_file
    df_ext = pd.read_csv(file_path, sep='\t', comment='!', index_col=0)

    with gzip.open(file_path, 'rt') as f:
        chars = None
        for line in f:
            if line.startswith('!Sample_characteristics_ch1'):
                chars = line.strip().split('\t')[1:]
                break

    # "Blood: Normal" -> 0, "Blood: CRC stage Duke A, B" / "Duke C, D" -> 1
    y_ext = np.array([0 if 'normal' in c.lower() else 1 for c in chars])
    print(f"GSE10715 — Samples: {len(y_ext)}, CRC: {sum(y_ext)}, Normal: {len(y_ext)-sum(y_ext)}")
    return df_ext, y_ext

X_ext, y_ext = load_gse10715()


# Align GSE10715 probes to the locked model's feature set


missing_probes = [p for p in final_selected_probes if p not in X_ext.index]
print(f"Probes missing in GSE10715: {len(missing_probes)} of {len(final_selected_probes)}")

X_ext_aligned = pd.DataFrame(index=final_selected_probes, columns=X_ext.columns, dtype=float)
for p in final_selected_probes:
    if p in X_ext.index:
        X_ext_aligned.loc[p] = X_ext.loc[p]
    else:

        X_ext_aligned.loc[p] = df.loc[p].mean()


#Predict and evaluate on GSE10715


probs_ext = final_model.predict_proba(X_ext_aligned.T)[:, 1]
preds_ext = (probs_ext >= 0.5).astype(int)

auc_ext = roc_auc_score(y_ext, probs_ext)
acc_ext = accuracy_score(y_ext, preds_ext)
sens_ext = recall_score(y_ext, preds_ext)
prec_ext = precision_score(y_ext, preds_ext, zero_division=0)
f1_ext = f1_score(y_ext, preds_ext, zero_division=0)
tn, fp, fn, tp = confusion_matrix(y_ext, preds_ext).ravel()
spec_ext = tn / (tn + fp)

print("\nGSE10715 TRUE EXTERNAL VALIDATION RESULTS")
print(f"AUC: {auc_ext:.3f}")
print(f"Accuracy: {acc_ext:.3f}")
print(f"Sensitivity: {sens_ext:.3f}")
print(f"Specificity: {spec_ext:.3f}")
print(f"Precision: {prec_ext:.3f}")
print(f"F1-score: {f1_ext:.3f}")

#Bootstrap 95% CI 
ci_ext_blood = bootstrap_ci_binary(y_ext, probs_ext, n_boot=2000)
print("\n 95% Bootstrap CI")
print(f"AUC 95% CI: [{ci_ext_blood['AUC_CI'][0]:.3f}, {ci_ext_blood['AUC_CI'][1]:.3f}]")
print(f"Accuracy 95% CI: [{ci_ext_blood['Accuracy_CI'][0]:.3f}, {ci_ext_blood['Accuracy_CI'][1]:.3f}]")
print(f"Sensitivity 95% CI: [{ci_ext_blood['Sensitivity_CI'][0]:.3f}, {ci_ext_blood['Sensitivity_CI'][1]:.3f}]")

# ROC figure
fpr_ext, tpr_ext, _ = roc_curve(y_ext, probs_ext)
plt.figure(figsize=(7,6))
plt.plot(fpr_ext, tpr_ext, color='green', lw=2, label=f"GSE10715 External (AUC={auc_ext:.3f})")
plt.plot([0,1],[0,1],'k--', alpha=0.4)
plt.xlabel("False Positive Rate"); plt.ylabel("True Positive Rate")
plt.title("External Validation on Independent Blood Cohort (GSE10715)")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "Figure_GSE10715_External_Validation_ROC.png"), dpi=300)
plt.show()

pd.DataFrame([{
    "Dataset": "GSE10715", "N": len(y_ext), "N_CRC": int(sum(y_ext)), "N_Normal": int(len(y_ext)-sum(y_ext)),
    "AUC": auc_ext, "Accuracy": acc_ext, "Sensitivity": sens_ext, "Specificity": spec_ext,
    "Precision": prec_ext, "F1": f1_ext
}]).to_csv(os.path.join(output_dir, "GSE10715_External_Validation_Results.csv"), index=False)








# BATCH-CORRECTION ATTEMPT: PER-STUDY Z-SCORE STANDARDIZATION
# uses only each dataset's OWN mean/SD, no labels involved,




# Re-train the final locked model using Z-SCORED GSE164191 values


# Compute per-probe mean/SD from GSE164191 
train_means = df.loc[final_selected_probes].mean(axis=1)
train_stds = df.loc[final_selected_probes].std(axis=1)

X_train_z = df.loc[final_selected_probes].sub(train_means, axis=0).div(train_stds, axis=0)

final_model_z = RandomForestClassifier(n_estimators=100, random_state=42)
final_model_z.fit(X_train_z.T, y)

print("Re-trained final model on z-scored GSE164191 values")


# Z-score GSE10715 using ITS OWN mean/SD 



ext_means = X_ext_aligned.mean(axis=1)
ext_stds = X_ext_aligned.std(axis=1)

X_ext_z = X_ext_aligned.sub(ext_means, axis=0).div(ext_stds, axis=0)


#Predict and evaluate again


probs_ext_bc = final_model_z.predict_proba(X_ext_z.T)[:, 1]
preds_ext_bc = (probs_ext_bc >= 0.5).astype(int)

auc_bc = roc_auc_score(y_ext, probs_ext_bc)
acc_bc = accuracy_score(y_ext, preds_ext_bc)
sens_bc = recall_score(y_ext, preds_ext_bc)
prec_bc = precision_score(y_ext, preds_ext_bc, zero_division=0)
f1_bc = f1_score(y_ext, preds_ext_bc, zero_division=0)
tn, fp, fn, tp = confusion_matrix(y_ext, preds_ext_bc).ravel()
spec_bc = tn / (tn + fp)

print("\n GSE10715 EXTERNAL VALIDATION — AFTER PER-STUDY Z-SCORE CORRECTION")
print(f"AUC: {auc_bc:.3f}")
print(f"Accuracy: {acc_bc:.3f}")
print(f"Sensitivity: {sens_bc:.3f}")
print(f"Specificity: {spec_bc:.3f}")
print(f"Precision: {prec_bc:.3f}")
print(f"F1-score: {f1_bc:.3f}")

ci_bc = bootstrap_ci_binary(y_ext, probs_ext_bc, n_boot=2000)
print("\n 95% Bootstrap CI (after correction)")
print(f"AUC 95% CI: [{ci_bc['AUC_CI'][0]:.3f}, {ci_bc['AUC_CI'][1]:.3f}]")
print(f"Accuracy 95% CI: [{ci_bc['Accuracy_CI'][0]:.3f}, {ci_bc['Accuracy_CI'][1]:.3f}]")
print(f"Sensitivity 95% CI: [{ci_bc['Sensitivity_CI'][0]:.3f}, {ci_bc['Sensitivity_CI'][1]:.3f}]")
print(f"Specificity 95% CI: not in your function — can add if needed")

# side-by-side comparison 
comparison_bc = pd.DataFrame({
    "Metric": ["AUC", "Accuracy", "Sensitivity", "Specificity", "Precision", "F1"],
    "Before_Correction": [auc_ext, acc_ext, sens_ext, spec_ext, prec_ext, f1_ext],
    "After_ZScore_Correction": [auc_bc, acc_bc, sens_bc, spec_bc, prec_bc, f1_bc]
})
print("\n", comparison_bc.to_string(index=False))
comparison_bc.to_csv(os.path.join(output_dir, "GSE10715_Batch_Correction_Comparison.csv"), index=False)
