from scipy.stats import mannwhitneyu
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import pandas as pd

#############################################################
# Input file, tab delimited
inf = 'example_in.txt'
df = pd.read_table(inf, index_col=0)

# How many categorial columns at the end of the table?
# This is important so that we don't try to test a categorial var
num_cat_columns = 2

# Which categorial variable is tested?
test_var = 'species'

# Which two categories are tested?
test_cat = ['red', 'blue']

#
outf = "MWU-test_corrected-pval.txt"
# Corrected p-value cutoff
# Anything above this will be ommited from output
cutoff = 0.05

which = "wilcoxon"  # or mwu
#############################################################



# collect results here
all_tests = []

# for each taxon make two lists, one for each category, amd test
for col in df.columns[:-1 * num_cat_columns]:
    catA = []
    catB = []
    for ind, row in df.iterrows():
        if row[test_var] == test_cat[0]:
            catA.append(row[col])
        elif row[test_var] == test_cat[1]:
            catB.append(row[col])
    s, p = None, None
    if which == "mwu":
        s, p = mannwhitneyu(catA, catB, use_continuity=True)
    elif which == "wilcoxon":
        s, p = wilcoxon(catA, y=catB, zero_method='wilcox', correction=False)

    all_tests.append([p, s, test_var, col])

# Add corrected p-values
tests = sorted(all_tests, reverse=True)

pvalues = [t[0] for t in tests]

stats = importr('stats')

p_adjust = stats.p_adjust(FloatVector(pvalues), method='BH')

tests_with_correct_pvals = zip(p_adjust, tests)

# Write output
headers = ['corrected P-values', 'p-values', 'statistic', 'factor', 'taxon']

outdf = pd.DataFrame([[i[0]] + i[1] for i in tests_with_correct_pvals if i[0] < cutoff],
                     columns=headers)

outdf.to_csv(outf, sep='\t')