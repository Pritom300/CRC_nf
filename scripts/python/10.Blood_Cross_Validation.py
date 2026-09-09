import sys
import pandas as pd
import numpy as np
import os
import gzip
from scipy import stats
from collections import Counter
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, auc, accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
import matplotlib.pyplot as plt
import warnings


if len(sys.argv) < 5:
    print("Error: Missing file arguments!")
    sys.exit(1)

importance_csv_file = sys.argv[1]
matrix_164191_file   = sys.argv[2]
gpl_platform_file    = sys.argv[3]
output_dir           = sys.argv[4]




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







warnings.filterwarnings("ignore")


def load_data():
    print("Loading GSE164191 (Blood Data)...")
    file_path = matrix_164191_file
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Dataset not found at {file_path}")

    df = pd.read_csv(file_path, sep='\t', comment='!', index_col=0)

    # Label extraction from titles
    with gzip.open(file_path, 'rt') as f:
        for line in f:
            if line.startswith('!Sample_title'):
                titles = line.strip().split('\t')[1:]
                break

    y = np.array([1 if any(x in t.lower() for x in ['crc', 'cancer', 'colorectal']) else 0 for t in titles])
    print(f"Samples: {len(y)}, Cancer: {sum(y)}, Healthy: {len(y)-sum(y)}")
    return df, y


# GENE MAPPING

def get_gene_mapping():
    file_path = gpl_platform_file
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"GPL570 mapping file not found at {file_path}")

    gpl = pd.read_csv(file_path, sep='\t', comment='#', low_memory=False)
    # Mapping Probe ID to a list of symbols
    mapping = {str(row['ID']): str(row['Gene Symbol']).split(' /// ') for _, row in gpl.iterrows()}
    return mapping


df = pd.read_csv(importance_csv_file)
target_genes = df["Gene"].tolist()


def select_features(X_train, y_train, mapping):
    selected = []
    cancer_idx = np.where(y_train == 1)[0]
    healthy_idx = np.where(y_train == 0)[0]

    for probe in X_train.index:
        symbols = mapping.get(probe, [])
        # Check if any symbol of this probe is in our target 72 genes
        if not any(s in target_genes for s in symbols):
            continue

        c = X_train.loc[probe].iloc[cancer_idx].values.astype(float)
        h = X_train.loc[probe].iloc[healthy_idx].values.astype(float)

        if len(c) < 2 or len(h) < 2:
            continue

        fc = np.mean(c) / (np.mean(h) + 1e-6)
        p = stats.ttest_ind(c, h, equal_var=False).pvalue

        # T-test and Fold Change filter
        if (fc > 1.1 or fc < 0.9) and p < 0.05:
            selected.append(probe)
    return selected


# BOOTSTRAP STABILITY

def bootstrap_feature_stability(X, y, mapping, n_boot=50, seed=123):
    rng = np.random.RandomState(seed)
    y = np.array(y)
    idx_cancer = np.where(y == 1)[0]
    idx_healthy = np.where(y == 0)[0]

    probe_hits = Counter()

    for b in range(n_boot):
        boot_idx = np.concatenate([
            rng.choice(idx_cancer, size=len(idx_cancer), replace=True),
            rng.choice(idx_healthy, size=len(idx_healthy), replace=True)
        ])
        X_boot = X.iloc[:, boot_idx]
        y_boot = y[boot_idx]

        selected = select_features(X_boot, y_boot, mapping)
        probe_hits.update(selected)

        if (b + 1) % 10 == 0:
            print(f"  bootstrap {b+1}/{n_boot} done")

    rows = []
    for p_id, freq in probe_hits.items():
        rows.append({
            "Probe_ID": p_id,
            "Gene_Symbol": ", ".join(mapping.get(p_id, ["N/A"])),
            "Times_Selected": freq,
            "Selection_Rate": freq / n_boot
        })

    freq_df = pd.DataFrame(rows).sort_values("Selection_Rate", ascending=False)
    freq_df.to_csv(os.path.join(output_dir, "Blood_Panel_Bootstrap_Stability.csv"), index=False)

    print(f"\nBLOOD-PANEL BOOTSTRAP STABILITY ({n_boot} resamples):")
    print(f"Mean selection rate: {freq_df['Selection_Rate'].mean()*100:.1f}%")
    print(f"Probes selected in >=80% of bootstraps: {(freq_df['Selection_Rate']>=0.8).sum()} / {len(freq_df)}")
    print(freq_df.head(15).to_string())

    return freq_df



def run_pipeline(X, y, mapping):
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    tprs, aucs, accs, precs, recs, f1s = [], [], [], [], [], []
    oof_true, oof_prob = [], []   # pooled out-of-fold storage
    all_selected_probes = []
    mean_fpr = np.linspace(0, 1, 100)

    plt.figure(figsize=(8, 6))

    for fold, (train_idx, test_idx) in enumerate(cv.split(X.T, y)):
        print(f"Processing Fold {fold+1}...")

        X_train, X_test = X.iloc[:, train_idx], X.iloc[:, test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Feature selection only on training data 
        selected = select_features(X_train, y_train, mapping)
        all_selected_probes.extend(selected)

        if not selected:
            continue

        model = RandomForestClassifier(n_estimators=100, random_state=42)
        model.fit(X_train.loc[selected].T, y_train)

        probs = model.predict_proba(X_test.loc[selected].T)[:, 1]
        oof_true.extend(y_test)      
        oof_prob.extend(probs)       
        fpr, tpr, _ = roc_curve(y_test, probs)

        tprs.append(np.interp(mean_fpr, fpr, tpr))
        aucs.append(auc(fpr, tpr))

        plt.plot(fpr, tpr, alpha=0.3, label=f'Fold {fold+1} (AUC = {auc(fpr, tpr):.2f})')

        preds = model.predict(X_test.loc[selected].T)
        accs.append(accuracy_score(y_test, preds))
        precs.append(precision_score(y_test, preds))
        recs.append(recall_score(y_test, preds))
        f1s.append(f1_score(y_test, preds))
        
        

    # Plotting Mean ROC 
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)

    plt.plot(mean_fpr, mean_tpr, color='blue', label=f'Mean ROC (AUC = {mean_auc:.3f} ± {std_auc:.3f})', lw=2)
    plt.plot([0, 1], [0, 1], linestyle='--', color='red', label='Random Guess')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Blood Validation ROC Curve (No Leakage)')
    plt.legend(loc="lower right")
    plt.grid(alpha=0.3)
    plt.savefig(os.path.join(output_dir, "Final_Blood_ROC.png"), dpi=300)
    plt.show()
    
    

    #  Metrics Output 
    
    print("\n" + "="*30)
    print("      FINAL RESULTS")
    print("="*30)
    print(f"AUC:       {np.mean(aucs):.3f} ± {std_auc:.3f}")
    print(f"Accuracy:  {np.mean(accs):.3f}")
    print(f"Precision: {np.mean(precs):.3f}")
    print(f"Recall:    {np.mean(recs):.3f}")
    print(f"F1-score:  {np.mean(f1s):.3f}")

    # pooled out-of-fold 95% CI 
    
    ci_blood = bootstrap_ci_binary(oof_true, oof_prob)

    print(f"AUC 95% CI:         [{ci_blood['AUC_CI'][0]:.3f}, {ci_blood['AUC_CI'][1]:.3f}]")
    print(f"Accuracy 95% CI:    [{ci_blood['Accuracy_CI'][0]:.3f}, {ci_blood['Accuracy_CI'][1]:.3f}]")
    print(f"Sensitivity 95% CI: [{ci_blood['Sensitivity_CI'][0]:.3f}, {ci_blood['Sensitivity_CI'][1]:.3f}]")
    print(f"Precision 95% CI:   [{ci_blood['Precision_CI'][0]:.3f}, {ci_blood['Precision_CI'][1]:.3f}]")
    print(f"F1 95% CI:          [{ci_blood['F1_CI'][0]:.3f}, {ci_blood['F1_CI'][1]:.3f}]")


    # Figure
    fpr8, tpr8, _ = roc_curve(oof_true, oof_prob)
    plt.figure(figsize=(7,6))
    plt.plot(fpr8, tpr8, color='darkorange', label=f"Pooled OOF (AUC={ci_blood['AUC_CI'].mean():.3f})")
    plt.plot([0,1],[0,1],'k--', alpha=0.4)
    plt.xlabel("False Positive Rate"); plt.ylabel("True Positive Rate")
    plt.title("Gene Refinement Validation (Blood, Leakage-Free 5-fold)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "Figure8_Blood_Refinement_ROC.png"), dpi=300)
    plt.show()


    # Figure End
    
    

    # Save All Genes to CSV 
    counts = Counter(all_selected_probes)
    gene_summary = []
    for p_id, freq in counts.items():
        gene_summary.append({
            "Probe_ID": p_id,
            "Gene_Symbol": ", ".join(mapping.get(p_id, ["N/A"])),
            "Selection_Frequency": freq
        })

    df_genes = pd.DataFrame(gene_summary).sort_values(by="Selection_Frequency", ascending=False)
    df_genes.to_csv(os.path.join(output_dir, "Validated_Genes_Blood.csv"), index=False)

    print("\nValidated gene list saved as 'Validated_Genes_Blood.csv'")
    print("\nTop 10 Most Consistent Genes in Blood:")
    print(df_genes.head(10))




if __name__ == "__main__":
    df, y = load_data()
    mapping = get_gene_mapping()
    run_pipeline(df, y, mapping)

    # bootstrap stability check
    blood_stability_df = bootstrap_feature_stability(df, y, mapping, n_boot=50)
