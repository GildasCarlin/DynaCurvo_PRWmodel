# ================================================================================================================
# ================================================================================================================
#                                                    packages importation
# ================================================================================================================
# ================================================================================================================

import numpy as np

import os
import regex as re
import pandas as pd

from scipy import stats
from scipy.stats import norm
import pingouin as pg
import scikit_posthocs as sp
import openturns as ot


import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from openturns.viewer import View

# ================================================================================================================
# ================================================================================================================
#                                                    intermediate functions
# ================================================================================================================
# ================================================================================================================

def checkNormality(quantity, getNorm=False, getBeta=False, getNPKS=False, showFig=True):
    # This function plot the Henry line of the given distribution
    # https://openturns.github.io/openturns/latest/auto_data_analysis/statistical_tests/plot_test_normality.html

    print('\n---------------------- check normality')

    histo_title = {'dx':r'$D_{\perp}$', 
                   'last_y':'Cell position at $T$=20 min (mm)',
                   'vd':'Velocity (µm/s)'}

    if quantity in ['last_y', 'r_d']:
        data_to_test = all_last.copy()
    elif quantity == 'dx':
        data_to_test = all_dx.copy()
    elif quantity == 'vd':
        data_to_test = all_vd.copy()
    #

    for assumption in assumptions:
        for celerity in celerities:
            print("\n...... assumption:", assumption, "celerity:", celerity, "µm/s")

            data = data_to_test[(data_to_test['assumption'] == assumption) & (data_to_test['celerity']==celerity)].copy()
            data = data[quantity]
        
            # compute normality test
            stat_ks = stats.kstest(data, "norm")
            print("... Kolmogorov–Smirnov: stat="+str(round(stat_ks.statistic, 3))+" pval="+str(stat_ks.pvalue)+" stat_sign="+str(stat_ks.statistic_sign))
            stat_and = stats.anderson(data, dist="norm")
            print("... Anderson-Darling: stat="+str(round(stat_and.statistic, 3))+"crit_val:", stat_and.critical_values)

            # open figure
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))
            txt = r"$\lambda$={}µm, $c$={}µm/s, ".format(int(wavelength), celerity)+ assumption + " cells, " + quantity
            plt.title(txt)

            # plot histogramm
            print("... plot histogramm")
            ax1.hist(data, density=True, color='blue', edgecolor='black')
        
            # plot pdf normal distribution
            x0, x1 = ax1.get_xlim()
            if getNorm:
                print("... plot normal distribution PDF")
    
                x = np.linspace(x0, x1, 100)
                pdf = norm.pdf(x, data.mean(), data.std())
                ax1.plot(x, pdf, color='red', linewidth=6.)
        
            # plot pdf beta distribution
            if getBeta:
                print("... plot Beta distribution PDF")
    
                pdf_graph = ot.Graph("", r"dx", "Density", True)
    
                bparameters = ot.BetaMuSigma(data.mean(), data.std(), 0., 1.)
                beta_dist = ot.ParametrizedDistribution(bparameters)
                pdf_beta = beta_dist.drawPDF().getDrawable(0)
                
                pdf_beta.setColor("#FEC6C5")
                pdf_beta.setLineWidth(5.)
                pdf_graph.add(pdf_beta)
                pdf_graph.setGrid(False)
                View(pdf_graph, figure=fig, axes=[ax1], add_legend=False)
            #

            # plot kernel smoothing pdf constructor 
            if getNPKS:
                print("... compute kernel smoothing methods for non-parametric distribution PDF")
                data_sample = ot.Sample(data.to_numpy().reshape(-1, 1))
                pdf_graph = ot.Graph("", r"dx", "Density", True)
    
                kernel = ot.KernelSmoothing(ot.Normal())
                bandwidth = kernel.computeMixedBandwidth(data_sample)
                kernel.setBoundaryCorrection(True)
                estimated = kernel.build(data_sample, bandwidth)
                kernel_plot = estimated.drawPDF().getDrawable(0)
    
                kernel_plot.setColor("#DA0176")
                kernel_plot.setLineWidth(5.)
                pdf_graph.add(kernel_plot)
                pdf_graph.setGrid(False)
                View(pdf_graph, figure=fig, axes=[ax1], add_legend=False)
            #

            # axes parameters
            ax1.set_xlim(left=x0, right=x1)
            ax1.set_xlabel(histo_title[quantity], fontsize=20)
            ax1.set_ylabel(r'Density', fontsize=20)
            ax1.tick_params(axis='both', which='major', labelsize=15)

            # legend
            handles, labels = plt.gca().get_legend_handles_labels()
            if quantity=='last_y' and celerity==.1:
                mpatch = mpatches.Patch(facecolor='red', label=r'Normal dist.', edgecolor='black')
                ax1.legend(handles=[mpatch], fontsize=15, loc="upper left")
            #

            if quantity=='dx' and celerity==.1:
                mpatch1 = mpatches.Patch(facecolor='red', label=r'Normal dist.', edgecolor='black')
                mpatch2 = mpatches.Patch(facecolor="#FEC6C5", label=r'Beta dist.', edgecolor='black')
                ax1.legend(handles=[mpatch1, mpatch2], fontsize=15, loc="upper left")
            #

            if quantity=='vd' and celerity==.1:
                mpatch = mpatches.Patch(facecolor='#DA0176', label=r'Non-parametric dist.', edgecolor='black')
                ax1.legend(handles=[mpatch], fontsize=15, loc="upper right")
            #


            # plot henry line
            pg.qqplot(data, dist='norm', ax=ax2, confidence=False)
            ax2.set_xlabel('Theoretical quantiles', fontsize=20)
            ax2.set_ylabel('Ordered quantiles', fontsize=20)
            ax2.tick_params(axis='both', which='major', labelsize=15)

            for text in ax2.texts:
                if re.match(r"\$R\^2=\d*.\d*\$", text.get_text()):
                    text.set_fontsize(18)
                #
            #    
        
            if showFig:
                plt.show()
            plt.close(fig)
    #
#

def checkCelerityEffect(quantity, assumption):
    ## select data
    factor = 'celerity'  ## independant variable
    dependant_variable = quantity 

    if dependant_variable in ['last_y', 'r_d']:
        data_to_test = all_last.copy()
    elif dependant_variable == 'dx':
        data_to_test = all_dx.copy()
    elif dependant_variable == 'vd':
        data_to_test = all_vd.copy()

    data_to_test = data_to_test[data_to_test['assumption']==assumption]
    print("---------- summary of data")
    print(data_to_test[[factor, dependant_variable]].groupby(factor).describe())
    #print(data_to_test[[factor, dependant_variable]].groupby(factor).median())
    #dx_df.groupby(['wavelength', 'celerity']).agg(["mean", "median", "std"]).reset_index()

    print('\n---------------------- Kruskal-Wallis test')
    # compute kruskal wallis test
    stat = pg.kruskal(
        data = data_to_test,
        dv = dependant_variable,
        between = factor,
    )
    print(stat)

    # compute  effect size 
    print('\n---------------------- Kruskal-Wallis complements')
    #Tomczak, E. & Tomczak, M. (2014). The need to report effect size estimates revisited. An overview of some recommended measures of effect size. TRENDS in Sport Sciences, 21(1).
    
    print('... eta square effect size through ANOVA')
    ranked_column = dependant_variable + '.rank'
    data_to_test[ranked_column] = data_to_test[dependant_variable].rank()
    stats_for_effect_size = pg.anova(data=data_to_test, dv=ranked_column, between=factor, effsize='n2')
    print('n^2 =', stats_for_effect_size.n2.values[0])

    print('... eta square effect size through formula')
    Hstatistic = stat.H.values[0]
    k = len(data_to_test[factor].unique())
    n = len(data_to_test[dependant_variable])
    
    eta2 = (Hstatistic-k+1.)/(n-k)
    print('n^2 =',eta2)

    print('... epsilon square effect size through formula')
    eps2 = Hstatistic*(n+1.)/(n**2-1.)
    print('eps^2 =',eps2)
    
    '''
    # compute ANOVA
    print('\n---------------------- ANOVA')
    print('... with partial n^2 effect size')
    aov1 = pg.anova(dv=dependant_variable, between=['assumption', 'celerity'], data=data_to_test, detailed=True)
    print(aov1)
    
    print('... with n^2 effect size')
    aov2 = pg.anova(dv=dependant_variable, between=['assumption', 'celerity'], data=data_to_test, detailed = True, effsize='n2')
    print(aov2)
    '''
    
    ### Dunn's test
    print("---------------------- post hoc Dunn's test")
    
    if dependant_variable in ['last_y', 'r_d']:
        data_to_test = all_last.copy()
    elif dependant_variable == 'dx':
        data_to_test = all_dx.copy()
    elif dependant_variable == 'vd':
        data_to_test = all_vd.copy()
    
    data_to_test = data_to_test[data_to_test['assumption']==assumption]
    # prepare data
    data_to_dunn = []
    for celerity in celerities:
        data_to_dunn.append(data_to_test[data_to_test[factor] == celerity][dependant_variable])
    
    p_values = sp.posthoc_dunn(data_to_dunn, p_adjust='holm')
    print(p_values)
    print(p_values <0.001)
    
    
    ### effect size
    print("\n---------------------- effect size")
    # Commun Language Effect Size (CLES): non paramétrique effect size https://doi.org/10.3758/s13428-015-0667-z
    all_d = []
    for data1 in data_to_dunn:
        line = []
        for data2 in data_to_dunn:
            eff=cles_safe(data1, data2)
            if eff<=.5:
                line.append(eff)
            else:
                line.append(1.-eff)
        all_d.append(line)
    
    print(pd.DataFrame(all_d))
#

def checkAssumptionEffect(quantity):
    factor = 'assumption'  ## independant variable
    dependant_variable = quantity

    if dependant_variable in ['last_y', 'r_d']:
        data_to_test = all_last.copy()
    elif dependant_variable == 'dx':
        data_to_test = all_dx.copy()
    elif dependant_variable == 'vd':
        data_to_test = all_vd.copy()

    #print(data_to_test[[factor, dependant_variable]].groupby(factor).mean())

    print('---------------------- Kruskal-Wallis test')
    # compute kruskal wallis test
    stat = pg.kruskal(
        data = data_to_test,
        dv = dependant_variable,
        between = factor,
    )
    print(stat)

    # compute  effect size 
    print('---------------------- Kruskal-Wallis complements')
    print('... eta square through formula')
    Hstatistic = stat.H.values[0]
    k = len(data_to_test[factor].unique())
    n = len(data_to_test[dependant_variable])
    
    eta2 = (Hstatistic-k+1.)/(n-k)
    print('n^2 =',eta2)
    
    print('... epsilon square through formula')
    eps2 = Hstatistic*(n+1.)/(n**2-1.)
    print('eps^2 =',eps2)
    # 

    ## compute welch t test
    print('\n---------------------- Welch t test between NDSC and DSC')

    for celerity in celerities:
        print('\n... celerity c={}µm/s'.format(celerity))
        data_to_welch=data_to_test[data_to_test['celerity']==celerity]
    
        data1 = data_to_welch[data_to_welch['assumption']=='NDSC'][dependant_variable]
        data2 = data_to_welch[data_to_welch['assumption']=='DSC'][dependant_variable]
        print("... taille d'échantillon,  NDSC : {} and DSC : {}".format(len(data1), len(data2)))
        print('... NDSC mean={}µm/s'.format(data1.mean()))
        print('... DSC mean={}µm/s'.format(data2.mean()))
        print('... Welch t test', stats.ttest_ind(data1, data2, equal_var=False))
        print("... Cohen's d =", pg.compute_effsize(data1, data2, eftype = 'cohen'))
        print('... Mann-Whitney U test', stats.mannwhitneyu(data1, data2))
        print("... CLM effect size =", pg.compute_effsize(data1, data2, eftype = 'CLES'))
    #
#

## different functions to compute CLES, validated and compared with pingouin package
def cles_mc(x, y, n=200_000):
    x = np.asarray(x)
    y = np.asarray(y)

    i = np.random.randint(0, len(x), n)
    j = np.random.randint(0, len(y), n)

    return np.mean(x[i] > y[j]) + 0.5 * np.mean(x[i] == y[j])
#

def cles_safe(x, y):
    x = np.asarray(x)
    y = np.asarray(y)

    # Tri pour approche efficace
    x = np.sort(x)
    y = np.sort(y)

    # Comptage efficace via recherche binaire
    idx = np.searchsorted(y, x, side='right')
    return idx.mean() / len(y)
#

# ================================================================================================================
# ================================================================================================================
#                                                    User interface
# ================================================================================================================
# ================================================================================================================

# study
substrate_regime="dyna"  # dyna or static
NDSC_study="NDSC_RightSinus_20min_VariousCelerity"
DSC_study="DSC_RightSinus_20min_VariousCelerity"
# assets
root=os.getcwd()
NDSC_assets=os.path.join(root, substrate_regime, "assets", NDSC_study)
DSC_assets=os.path.join(root, substrate_regime, "assets", DSC_study)

assumptions=['NDSC', 'DSC']
mapping_path_assumption={"NDSC":NDSC_assets, "DSC":DSC_assets}

## conditions experimentales
wavelength = 20.

## list of wave celerity to plot distributions
celerities = [0., .1, .7, .9, 5.]

# what to do
normality=False
celerity_effect=True
assumption_effect=False

# what to test
quantity="last_y"
nature="NDSC"

# ================================================================================================================
# ================================================================================================================
#                                                    MAIN function
# ================================================================================================================
# ================================================================================================================

print('\n------------------------------------ get general indicators')
#print('---------------------- read experimental data')
#dx_exp = pd.read_csv(root + "\\experimental_static\\summary_dx_sin20.csv")
#vd_exp = pd.read_csv(root + "\\experimental_static\\summary_vd_sin20.csv")

print('---------------------- read numerical data')
#### NS
assumption = 'NDSC'
last_NS = pd.read_csv(NDSC_assets + "\\" + "last.csv")
last_NS = last_NS[last_NS['celerity'].isin(celerities)]
last_NS = last_NS[last_NS['wavelength']==wavelength]
last_NS.insert(0, 'assumption', len(last_NS) * [assumption])

dx_NS = pd.read_csv(NDSC_assets + "\\" + "dx.csv")
dx_NS = dx_NS[dx_NS['celerity'].isin(celerities)]
dx_NS = dx_NS[dx_NS['wavelength']==wavelength]
dx_NS.insert(0, 'assumption', len(dx_NS) * [assumption])

vd_NS = pd.read_csv(NDSC_assets + "\\" + "vd.csv")
vd_NS = vd_NS[vd_NS['celerity'].isin(celerities)]
vd_NS = vd_NS[vd_NS['wavelength']==wavelength]
vd_NS.insert(0, 'assumption', len(vd_NS) * [assumption])

#### S
assumption = 'DSC'
last_S = pd.read_csv(DSC_assets + "\\" + "last.csv")
last_S = last_S[last_S['celerity'].isin(celerities)]
last_S = last_S[last_S['wavelength']==wavelength]
last_S.insert(0, 'assumption', len(last_S) * [assumption])

dx_S = pd.read_csv(DSC_assets + "\\" + "dx.csv")
dx_S = dx_S[dx_S['celerity'].isin(celerities)]
dx_S = dx_S[dx_S['wavelength']==wavelength]
dx_S.insert(0, 'assumption', len(dx_S) * [assumption])

vd_S = pd.read_csv(DSC_assets + "\\" + "vd.csv")
vd_S = vd_S[vd_S['celerity'].isin(celerities)]
vd_S = vd_S[vd_S['wavelength']==wavelength]
vd_S.insert(0, 'assumption', len(vd_S) * [assumption])

######## concat
all_last = pd.concat([last_NS, last_S], ignore_index=True)
all_dx = pd.concat([dx_NS, dx_S], ignore_index=True)
all_vd = pd.concat([vd_NS, vd_S], ignore_index=True)

### convert to mm
all_last['last_y'] /= 1000.

if normality:
    checkNormality(quantity=quantity, getNorm=True)

if celerity_effect:
    checkCelerityEffect(quantity=quantity, assumption=nature)

if assumption_effect:
    checkAssumptionEffect(quantity=quantity)






"""
# ================================================================================================================
#                                                Comparison with static exp distribution
# ================================================================================================================
## select data
factor = 'assumption'  ## independant variable
dependant_variable = 'dx'

if dependant_variable in ['last_y', 'r_d']:
    data_to_test = all_last.copy()
elif dependant_variable == 'dx':
    data_to_test = all_dx.copy()
    data_exp = dx_exp[dependant_variable].copy()
elif dependant_variable == 'vd':
    data_to_test = all_vd.copy()
    data_exp = vd_exp[dependant_variable].copy()


print('\n---------------------- Welch t test between static exp and DSC / NDSC')

for celerity in celerities:
    print('\n... celerity c={}µm/s'.format(celerity))
    data_to_welch=data_to_test[data_to_test['celerity']==celerity]

    data1 = data_to_welch[data_to_welch['assumption']=='NDSC'][dependant_variable]
    data2 = data_to_welch[data_to_welch['assumption']=='DSC'][dependant_variable]
    print("... taille d'échantillon,  NDSC: {}, DSC: {} and T-cells: {}".format(len(data1), len(data2), len(data_exp)))
    print('... NDSC mean={}µm/s'.format(data1.mean()))
    print('... DSC mean={}µm/s'.format(data2.mean()))
    print('... T-cells mean={}µm/s'.format(data_exp.mean()))

    print('... NDSC vs exp:', stats.ttest_ind(data1, data_exp, equal_var=False))
    print("... Cohen's d =", pg.compute_effsize(data1, data_exp, eftype = 'cohen'))
    
    print('... DSC vs exp:', stats.ttest_ind(data2, data_exp, equal_var=False))
    print("... Cohen's d =", pg.compute_effsize(data2, data_exp, eftype = 'cohen'))
#

"""


