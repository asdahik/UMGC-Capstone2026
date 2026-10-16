import pandas as pd, numpy as np, csv
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
from sklearn.inspection import permutation_importance
from imputation import PCAImputer as imp


# classification labels
LEVELS = ["low", "intermediate", "high"]
# ------- Load ------

# get the working/current directory
wd = Path.cwd()/"Metabolomics"
with open(wd / "cultivar_metabolome.txt", encoding="UTF-8", errors="replace") as f:
    rows = list(csv.reader(f, delimiter="\t"))
header, data = rows[12], rows[13:]
mst = header[10:] # 115 metabolite columns

df = pd.DataFrame(data, columns=header)
df[mst] = df[mst].replace("NA", np.nan).apply(pd.to_numeric)
df["Cultivar ID"] = df["Cultivar ID"].astype(int)

# annotation lookup: MST code -> readable name
names = dict(zip(mst, rows[3][10:]))
label = {m: (names[m] if names[m] != "NA" else m) for m in mst}

tolerance = pd.read_csv(f"{wd}/cultivar_drought_tolerance.txt", sep="\t", skiprows=1).iloc[:, :4]
tolerance.columns = ["Cultivar ID", "Cultivar", "breeder", "DRYM"] # ignore breeder column
tolerance["tol"] = pd.qcut(tolerance.DRYM, 3, labels=LEVELS)
df = df.merge(tolerance[["Cultivar ID", "DRYM", "tol"]], on="Cultivar ID")

# train random forest
train = df["Trial ID"].str.startswith("F")
X_tr, X_va = df.loc[train, mst], df.loc[-train, mst]
y_tr = df.loc[train, "tol"].astype(str).values
y_va = df.loc[~train, "tol"].astype(str).values
g_tr = df.loc[train, "Cultivar ID"].values
g_va = df.loc[~train, "Cultivar ID"].values

# impute (fit on training only), then:
# rf = RandomForestClassifier(n_estimators=1000, oob_score=True, n_jobs=-1, random_state=1)
# rf.fit(X_train, df.loc[train, "tol"])

# ------- Run Imputation -------
imputation = imp(k=5).fit(X_tr.values)
Xtr = imputation.train_filled_
Xva = imputation.transform(X_va.values)

# ------- Train + validate

rf = RandomForestClassifier(n_estimators=1000, oob_score=True, n_jobs=-1, random_state=1)
rf.fit(Xtr, y_tr)

print(f"OOB error (training trials):\t\t{1-rf.oob_score_:.3f}")
print(f"Accuracy on agronomic validation:\t\t{(rf.predict(Xva) == y_va).mean():.3f}")
print(pd.DataFrame(confusion_matrix(y_va, rf.predict(Xva), labels=LEVELS),index=[f"obs_{l}" for l in LEVELS], columns=LEVELS))

# with the results being outputted from RandomForestClassifier, write output to a results file
output = open(f"{wd}/metabolome_output.txt", "w")
output.write("Metabolite\tgini_coefficient")
gini = pd.Series(rf.feature_importances_, index=[label[m] + f" ({m})" for m in mst])
print("\nTop 15 metabolites (Gini importance):")
print(gini.sort_values(ascending=False).head(15).round(4))
output.write((gini.sort_values(ascending=False).head(15).round(4)).to_string())
# ------ Honest estimate: hold out whole cultivars ------
# Each cultivar appears in bot train and validation with the same label, so 
# the numbers above should partly reflect recognising the cultivar. Held-out
# cultivars are absent from training (imputer refit inside each fold)

all_cv = np.unique(np.r_[g_tr, g_va])
folds = np.array_split(np.random.RandomState(0).permutation(all_cv), 5)

# prediction, observation
pred, obs = [], []
for held in folds:
    tr = ~np.isin(g_tr, held)
    va = np.isin(g_va, held)
    im = imp(k=5).fit(X_tr.values[tr])
    m = RandomForestClassifier(n_estimators=500, n_jobs=-1, random_state=1)
    m.fit(im.train_filled_, y_tr[tr])
    pred += list(m.predict(im.transform(X_va.values[va])))
    obs += list(y_va[va])

pred, obs = np.array(pred), np.array(obs)
base = pd.Series(y_va).value_counts(normalize=True).max()
print(f"\nLeave-cultivar-out accuracy (agronomic): {(pred == obs).mean():.3f}"
      f"   | majority-class baseline: {base:.3f}")
output.write(f"\nLeave-cultivar-out accuracy (agronomic): {(pred == obs).mean():.3f}"
      f"   | majority-class baseline: {base:.3f}")

output.close()
# output results 
# Optional: permutation importance on held-out validation data (less biased than Gini)
# pi = permutation_importance(rf, Xva, y_va, n_repeats=10, n_jobs=-1, random_state=1)