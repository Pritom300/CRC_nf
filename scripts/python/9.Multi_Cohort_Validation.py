import sys
import os
import pandas as pd
import numpy as np
import re
import gzip
import json
import warnings
warnings.filterwarnings("ignore")

from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve, accuracy_score, recall_score, precision_score, f1_score, confusion_matrix
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import shap


if len(sys.argv) < 9:
    print("Error: Missing file arguments!")
    sys.exit(1)

importance_csv_file = sys.argv[1]
train_matrix_file    = sys.argv[2]
matrix_110223_file   = sys.argv[3]
matrix_164191_file   = sys.argv[4]
matrix_142987_file   = sys.argv[5]
gpl_platform_file    = sys.argv[6]
ensg_json_file       = sys.argv[7]
output_dir           = sys.argv[8]




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





# TARGET GENES


csv_path = importance_csv_file

if not os.path.exists(csv_path):
    raise FileNotFoundError(f"CSV did not find: {csv_path}")


df = pd.read_csv(csv_path)


if "Gene" not in df.columns:
    raise ValueError(
        f"CSV-has not gene Columns: {list(df.columns)}"
    )


# Duplicate will remove
target_genes = list(dict.fromkeys(
    df["Gene"]
    .dropna()
    .astype(str)
    .str.strip()
    .str.upper()
))

print(f"Total target gene: {len(target_genes)}")
print(target_genes)



def load_gse41258():
    file_path = train_matrix_file
    titles, chars = [], []

    with gzip.open(file_path, 'rt') as f:
        for line in f:
            if line.startswith('!Sample_title'):
                titles = [t.strip('"') for t in line.split('\t')[1:]]
            elif line.startswith('!Sample_characteristics_ch1'):
                chars.append([t.strip('"') for t in line.split('\t')[1:]])

    df = pd.read_csv(file_path, sep='\t', comment='!', index_col=0)

    chars_df = pd.DataFrame(chars).T.fillna('')
    combined = chars_df.apply(lambda x: ' '.join(x).lower(), axis=1).tolist()
    texts = [f"{t.lower()} {c}" for t, c in zip(titles, combined)]

    labels = []
    for txt in texts:
        if re.search(r"(metastasis|primary|carcinoma|crc)", txt):
            labels.append(1)
        elif re.search(r"(normal|control|healthy)", txt):
            labels.append(0)
        else:
            labels.append(-1)

    idx = [i for i,x in enumerate(labels) if x!=-1]
    return df.iloc[:,idx], [labels[i] for i in idx]
    
    

def load_gse110223():
    titles=[]
    with gzip.open(matrix_110223_file,'rt') as f:
        for line in f:
            if line.startswith('!Sample_title'):
                titles=[t.strip('"') for t in line.split('\t')[1:]]
                break
    df=pd.read_csv(matrix_110223_file,sep='\t',comment='!',index_col=0)
    labels=[1 if 'cancer' in t.lower() else 0 for t in titles]
    return df,labels
    

def load_gse164191():
    df=pd.read_csv(matrix_164191_file,sep='\t',comment='!',index_col=0)
    return df,[0]*62+[1]*59
    

def load_gse142987():
    df=pd.read_csv(matrix_142987_file,sep='\t',index_col=0)
    labels=[1 if c.startswith('L') else 0 for c in df.columns]
    return df,labels
    

def load_mappings():
    gpl=pd.read_csv(gpl_platform_file,sep='\t',comment='#',low_memory=False)
    gpl_map={r['ID']:str(r['Gene Symbol']).split(' /// ')[0] for _,r in gpl.iterrows()}
    with open(ensg_json_file) as f:
        ensg=json.load(f)
    return gpl_map,{v:k for k,v in ensg.items()}


def preprocess(df, mapping, scaler=None, fit=False, is_ensg=False):
    if is_ensg:
        df.index=[i.split('.')[0] for i in df.index]

    df['Symbol']=df.index.map(mapping)
    df=df.dropna(subset=['Symbol'])

    X=df.groupby('Symbol').mean().T

    for g in target_genes:
        if g not in X.columns:
            X[g]=0
    X=X[target_genes]

    if X.max().max()>50:
        X=np.log2(X+1)

    if fit:
        scaler=StandardScaler()
        X=scaler.fit_transform(X)
    else:
        X=scaler.transform(X)

    return pd.DataFrame(X,columns=target_genes,index=X.index if isinstance(X,pd.DataFrame) else None), scaler


# METRICS

def evaluate(y, probs, name):
    pred=(probs>=0.5).astype(int)

    tn,fp,fn,tp=confusion_matrix(y,pred).ravel()

    print(f"\n{name}")
    print("AUC:",roc_auc_score(y,probs))
    print("Accuracy:",accuracy_score(y,pred))
    print("Sensitivity:",tp/(tp+fn))
    print("Specificity:",tn/(tn+fp))
    print("F1:",f1_score(y,pred))
    
    # bootstrap 95% CI 
    ci_result = bootstrap_ci_binary(y, probs)
    print(f"  AUC 95% CI:         [{ci_result['AUC_CI'][0]:.3f}, {ci_result['AUC_CI'][1]:.3f}]")
    print(f"  Accuracy 95% CI:    [{ci_result['Accuracy_CI'][0]:.3f}, {ci_result['Accuracy_CI'][1]:.3f}]")
    print(f"  Sensitivity 95% CI: [{ci_result['Sensitivity_CI'][0]:.3f}, {ci_result['Sensitivity_CI'][1]:.3f}]")
    print(f"  Precision 95% CI:   [{ci_result['Precision_CI'][0]:.3f}, {ci_result['Precision_CI'][1]:.3f}]")
    print(f"  F1 95% CI:          [{ci_result['F1_CI'][0]:.3f}, {ci_result['F1_CI'][1]:.3f}]")


gpl_map,ensg_map=load_mappings()

# TRAIN
df_train,y_train=load_gse41258()
X_train,scaler=preprocess(df_train,gpl_map,fit=True)

model=RandomForestClassifier(n_estimators=100,random_state=42)
model.fit(X_train,y_train)

# TEST
datasets=[]
for name,loader,args in [
    ("Tissue",load_gse110223,(gpl_map,False)),
    ("Blood",load_gse164191,(gpl_map,False)),
    ("Plasma",load_gse142987,(ensg_map,True))
]:
    df,y=loader()
    X,_=preprocess(df,args[0],scaler=scaler,is_ensg=args[1])
    p=model.predict_proba(X)[:,1]
    evaluate(y,p,name)
    datasets.append((name,y,p,X))

# Figure
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve

plt.figure(figsize=(7,6))
for name, y, p, _ in datasets:   
    fpr, tpr, _ = roc_curve(y, p)
    auc_val = roc_auc_score(y, p)
    plt.plot(fpr, tpr, label=f"{name} (AUC={auc_val:.3f})")
plt.plot([0,1],[0,1],'k--', alpha=0.4)
plt.xlabel("False Positive Rate"); plt.ylabel("True Positive Rate")
plt.title("Validation in Tissue, Blood, and Plasma")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "Figure7_Tissue_Blood_Plasma_ROC.png"), dpi=300)
plt.show()

# Figure End


# ROC PLOT

plt.figure(figsize=(7,6))
for name,y,p,_ in datasets:
    fpr,tpr,_=roc_curve(y,p)
    plt.plot(fpr,tpr,label=name)
plt.plot([0,1],[0,1],'k--')
plt.legend()
plt.title("ROC Comparison")
plt.savefig(os.path.join(output_dir, "ROC.png"),dpi=300)
plt.show()


# SHAP Analysis

explainer=shap.TreeExplainer(model)
shap_vals = explainer.shap_values(X_train)

# handle all SHAP formats safely
if isinstance(shap_vals, list):
    shap_vals = shap_vals[1]   # class 1
elif len(shap_vals.shape) == 3:
    shap_vals = shap_vals[:, :, 1]  # class 1

# now shape = (samples, features)
mean_shap = np.abs(shap_vals).mean(axis=0)

top_genes = pd.Series(mean_shap, index=X_train.columns)\
    .sort_values(ascending=False)\
    .head(10)

print("\nTop Genes:",top_genes.index.tolist())

plt.figure(figsize=(8,5))
top_genes.plot(kind='bar')
plt.title("Top 10 SHAP Genes")
plt.savefig(os.path.join(output_dir, "SHAP_top10.png"),dpi=300)
plt.show()
