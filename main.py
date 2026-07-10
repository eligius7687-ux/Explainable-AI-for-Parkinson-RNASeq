import pandas as pd
import numpy as np

from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests


df = pd.read_csv(
    "data/GSE165082_PD-CC.counts.txt",
    sep="\t"
)

# جدا کردن داده‌ها
genes = df["Geneid"]

counts = df.drop("Geneid", axis=1)

# Filtering
filtered_counts = counts[counts.sum(axis=1) >= 10]

# Normalization
normalized_counts = np.log2(filtered_counts + 1)


# گروه‌ها
pd_samples = [col for col in normalized_counts.columns if "PD" in col]
cc_samples = [col for col in normalized_counts.columns if "CC" in col]


# آزمون t-test برای هر ژن

p_values = []

for gene in normalized_counts.index:
    
    pd_values = normalized_counts.loc[gene, pd_samples]
    cc_values = normalized_counts.loc[gene, cc_samples]
    
    stat, p = ttest_ind(
        pd_values,
        cc_values,
        equal_var=False
    )
    
    p_values.append(p)


# ساخت جدول نتایج

results = pd.DataFrame({
    "Gene": genes[filtered_counts.index],
    "p_value": p_values
})


# اصلاح FDR
results["FDR"] = multipletests(
    results["p_value"],
    method="fdr_bh"
)[1]


# مرتب‌سازی
results = results.sort_values(
    "FDR"
)


print(results.head(20))