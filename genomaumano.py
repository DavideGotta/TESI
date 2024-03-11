import GEOparse as geo
import scipy as sc
import scipy.stats
import numpy as np

df= geo.get_GEO(filepath='GSE4040_family.soft.gz')
gse=0
gpls=0
#ns, ng = df.shape
print(df.gsms, df.gpls)


df.iloc[[1, 2, 3], ng - 1] = 'Wild type'
df.iloc[[4, 5, 6], ng - 1] = 'Inhibitor'
print(df)

G1 = df.loc[df.CLASS == 'Wild type', df.columns != 'CLASS']
print(G1)
G2 = df.loc[df.CLASS == 'Inhibitor', df.columns != 'CLASS']
print(G2)

st, pv = sc.stats.ttest_ind(G1, G2)
pvc = sc.stats.false_discovery_control(pv)

alpha = 0.05
col_test = G1.columns[pv < alpha]
gene_names_ttest_gsea = geo.gene_names(col_test, gpls, 'Gene Symbol')
print('T-test results for GSEA,(' + str(len(gene_names_ttest_gsea)) + '):')
for gene_ttest_gsea in gene_names_ttest_gsea:
    print(gene_ttest_gsea)

alpha = 0.05
col_test = G1.columns[pv < alpha]
gene_names_ttest_manuale = geo.gene_names(col_test, gpls, 'Gene Symbol')
print('T-test results (manuale), (' + str(len(gene_names_ttest_manuale)) + '):')
for gene_ttest_manuale in gene_names_ttest_manuale:
    print(gene_ttest_manuale)

df.iloc[:, df.columns != 'CLASS'] = 2 ** (df.iloc[:, df.columns != 'CLASS'])

means = df.groupby('CLASS').mean()
dm = means.iloc[0] / means.iloc[1]
log10pv = -1 * np.log10(pv)
log2fc = np.log2(dm)
down, up = geo.volcano_plot(log2fc, log10pv, thres=[-2, 2, 1.8])
gene_names_down_gsea = geo.gene_names(down, gpls, 'Gene Symbol')
print('Down Volcano results for GSEA, (' + str(len(gene_names_down_gsea)) + '):')
for gene_down_gsea in gene_names_down_gsea:
    print(gene_down_gsea)

gene_names_up_gsea = geo.gene_names(up, gpls, 'Gene Symbol')
print('Up Volcano results for GSEA, (' + str(len(gene_names_up_gsea)) + '):')
for gene_up_gsea in gene_names_up_gsea:
    print(gene_up_gsea)

down, up = geo.volcano_plot(log2fc, log10pv, thres=[-3, 3.5, 3.5])
gene_names_down_manuale = geo.gene_names(down, gpls, 'Gene Symbol')
print('Down Volcano results (manuale), (' + str(len(gene_names_down_manuale)) + '):')
for gene_down_manuale in gene_names_down_manuale:
    print(gene_down_manuale)
gene_names_up_manuale = geo.gene_names(up, gpls, 'Gene Symbol')
print('Up Volcano results (manuale), (' + str(len(gene_names_up_manuale)) + '):')
for gene_up_manuale in gene_names_up_manuale:
    print(gene_up_manuale)

st1, pv1 = scipy.stats.mannwhitneyu(G1, G2)
print(st1)
print(pv1)
col_utest = G1.columns[pv1 < 2.2 * 10 ** (-3)]
print(col_utest)
gene_names_utest_gsea = geo.gene_names(col_utest, gpls, 'Gene Symbol')

print('U-test results for GSEA, (' + str(len(gene_names_utest_gsea)) + '):')
for gene_utest_gsea in gene_names_utest_gsea:
    print(gene_utest_gsea)

st1, pv1 = scipy.stats.mannwhitneyu(G1, G2)
print(st1)
print(pv1)
col_utest = G1.columns[pv1 < 2.2 * 10 ** (-3)]
print(col_utest)
gene_names_utest_manuale = geo.gene_names(col_utest, gpls, 'Gene Symbol')

print('U-test results (manuale), (' + str(len(gene_names_utest_manuale)) + '):')
for gene_utest_manuale in gene_names_utest_manuale:
    print(gene_utest_manuale)

st2, pv2 = scipy.stats.f_oneway(G1, G2)
print(st2)
print(pv2)
col_anova = G1.columns[pv2 < 0.005]
gene_names_anova_gsea = geo.gene_names(col_anova, gpls, 'Gene Symbol')
print('ANOVA results for GSEA, (' + str(len(gene_names_anova_gsea)) + '):')
for gene_anova_gsea in gene_names_anova_gsea:
    print(gene_anova_gsea)

st2, pv2 = scipy.stats.f_oneway(G1, G2)
print(st2)
print(pv2)
col_anova = G1.columns[pv2 < 0.000005]
gene_names_anova_manuale = geo.gene_names(col_anova, gpls, 'Gene Symbol')
print('ANOVA results (manuale),(' + str(len(gene_names_anova_manuale)) + '):')
for gene_anova_manuale in gene_names_anova_manuale:
    print(gene_anova_manuale)

st3, pv3 = scipy.stats.kruskal(G1, G2)
print(st3)
print(pv3)
col_kruskal = G1.columns[pv3 < 0.009]
print(col_kruskal)
gene_names_kruskal_gsea = geo.gene_names(col_kruskal, gpls, 'Gene Symbol')
print('KRUSKAL results for GSEA, (' + str(len(gene_names_kruskal_gsea)) + '):')
for gene_kruskal_gsea in gene_names_kruskal_gsea:
    print(gene_kruskal_gsea)

st3, pv3 = scipy.stats.kruskal(G1, G2)
print(st3)
print(pv3)
col_kruskal = G1.columns[pv3 < 0.004]
print(col_kruskal)
gene_names_kruskal_manuale = geo.gene_names(col_kruskal, gpls, 'Gene Symbol')
print('KRUSKAL results (manuale), (' + str(len(gene_names_kruskal_manuale)) + '):')
for gene_kruskal_manuale in gene_names_kruskal_manuale:
    print(gene_kruskal_manuale)