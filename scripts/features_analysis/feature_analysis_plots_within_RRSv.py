# This script automates the exploratory analysis of the PRS and the individual replicates of a version of RRS.
# Author: Chop Yan Lee

import sys
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

all_features= ['Probability', 'IUPredShort', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLC', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'DomainFreqbyProtein1', 'DomainFreqinProteome1']
all_features_renamed= ['Probability', 'IUPredShort', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLC', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'DomainFreqbyProtein', 'DomainFreqinProteome']
all_features_reordered= ['Probability', 'IUPredShort', 'Anchor', 'DomainOverlap', 'DomainEnrichment_pvalue', 'qfo_RLC', 'vertebrates_RLC', 'mammalia_RLC',  'metazoa_RLC', 'qfo_RLCvar', 'vertebrates_RLCvar', 'mammalia_RLCvar', 'metazoa_RLCvar', 'DomainEnrichment_zscore', 'DomainFreqbyProtein', 'DomainFreqinProteome']
all_features_reordered_label= ['Probability', 'IUPredShort', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'vertebrates_RLC', 'mammalia_RLC',  'metazoa_RLC', 'qfo_RLCvar', 'vertebrates_RLCvar', 'mammalia_RLCvar', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'DomainFreqbyProtein', 'DomainFreqinProteome', 'label']

DMI_count_df= pd.DataFrame(data= {'Class': ['CLV', 'DEG', 'DOC', 'LIG', 'MOD', 'TRG'], 'ElmDB': [11, 25, 31, 165, 37, 22]})

fontsize= 12
title_fontsize= 14
sns.set_style('darkgrid')

def preprocessing_dataset(PRS_input, RRS_input_list): 
    """
    Preprocess the NaNs and dummy values in the PRS and RRS replicates. Rows with NaN are removed. Dummy value (88888) is replaced with the median of the feature in the dataset, i.e. median of the feature in PRS or median in RRS. All RRS replicates are concatenated into one single RRS dataframe.
    
    Args:
        PRS_path (str): Absolute path to the PRS dataset
        RRS_input_list (list of str): List of absolute path to the individual replicates of an RRS version

    Returns:
        PRS (pd.DataFrame): The processed PRS
        RRS (pd.DataFrame): The processed RRS (all triplicates concatenated into one RRS dataframe)
    """
    PRS= pd.read_csv(PRS_input, sep= '\t', index_col= 0)
    PRS.replace(88888, PRS.DomainEnrichment_zscore.median(), inplace= True)
    RRS= pd.DataFrame(columns= PRS.columns)
    for RRS_input in RRS_input_list:
        df= pd.read_csv(RRS_input, sep= '\t', index_col= 0)
        df.replace(88888, df.DomainEnrichment_zscore.median(), inplace= True)
        RRS= pd.concat([RRS, df], axis= 0)
    for df in [PRS, RRS]:
        for ind, row in df.iterrows():
            if pd.notna(row['DomainFreqbyProtein2']):
                df.loc[ind, 'DomainFreqbyProtein1'] = np.mean([row['DomainFreqbyProtein1'], row['DomainFreqbyProtein2']])
                df.loc[ind, 'DomainFreqinProteome1'] = np.mean([row['DomainFreqinProteome1'], row['DomainFreqinProteome2']])
        df.dropna(subset= all_features, inplace= True)
        df.rename(columns= {'DomainFreqbyProtein1': 'DomainFreqbyProtein', 'DomainFreqinProteome1': 'DomainFreqinProteome'}, inplace= True)
    return PRS, RRS

def make_DMI_fraction_plot(PRS, RRS):
    """
    Plot the representation of the DMI types in each ELM classes (DEG, DOC, etc.) in the PRS and RRS using the fraction of DMI type in each ELM classes. The plot is saved as pdf file.

    Args:
        PRS (pd.DataFrame): The processed PRS from the function preprocessing_dataset
        RRS (pd.DataFrame): The processed RRS (all triplicates concatenated into one RRS dataframe) from the function preprocessing_dataset
    """
    global DMI_count_df
    for i, df in enumerate([PRS, RRS]):
        if i == 0:
            temp= pd.DataFrame(data= {'Class': pd.Series(df.Elm.unique()).str.slice(stop= 3).value_counts().index, 'PRS': pd.Series(df.Elm.unique()).str.slice(stop= 3).value_counts().values})
        else:
            temp= pd.DataFrame(data= {'Class': pd.Series(df.Elm.unique()).str.slice(stop= 3).value_counts().index, 'RRS': pd.Series(df.Elm.unique()).str.slice(stop= 3).value_counts().values})
        DMI_count_df= DMI_count_df.merge(temp, how= 'inner')
    DMI_count_df['PRS_fraction']= DMI_count_df['PRS']/DMI_count_df['ElmDB']
    DMI_count_df['RRS_fraction']= DMI_count_df['RRS']/DMI_count_df['ElmDB']

    N= 6
    ind= np.arange(N)
    width= 0.3

    plt.figure(figsize= (8,6))
    plt.bar(ind, DMI_count_df.PRS_fraction, width, color= 'c', label= 'PRS')
    plt.bar(ind + width, DMI_count_df.RRS_fraction, width, color= 'b', label= 'RRS')

    plt.xticks(ind + width / 2, DMI_count_df.Class)
    plt.title(f'Fraction of DMI represented in the PRS and {RRS_version} by class', fontsize= title_fontsize)
    plt.ylabel('Fraction of DMI types with DMI instance', fontsize= fontsize)
    plt.legend(loc= 'best')
    # plt.grid(alpha= 0.2)
    plt.ylim([0, 1.0])
    plt.savefig(plot_path + f'/DMI_fraction_PRS_{RRS_version}_1_2_3.pdf', bbox_inches= 'tight')
    print(f'DMI fraction plot of {RRS_version} saved in {plot_path}.')
    plt.close()

def make_feature_violin_plots(PRS, RRS):
    """
    Plot the distribution of different features in the PRS and RRS using violinplots. The plots are saved as pdf file.

    Args:
        PRS (pd.DataFrame): The processed PRS from the function preprocessing_dataset
        RRS (pd.DataFrame): The processed RRS (all triplicates concatenated into one RRS dataframe) from the function preprocessing_dataset
    """
    fig, ax= plt.subplots(figsize= (10, 6))

    PRS_data= [PRS[feature] for feature in all_features_reordered[1:-3]]
    RRS_data= [RRS[feature] for feature in all_features_reordered[1:-3]]

    vp1= ax.violinplot(PRS_data, showmedians= True, widths= 0.8)
    for pc in vp1['bodies']:
        pc.set_edgecolor('black')
    vp2= ax.violinplot(RRS_data, showmedians= True, widths= 0.8)

    plt.setp(ax, xticks= [x + 1 for x in range(len(all_features_reordered[1:-3]))], xticklabels= all_features_reordered[1:-3])
    plt.setp(ax.xaxis.get_majorticklabels(), rotation = 45, horizontalalignment= 'right')
    plt.legend([vp1['bodies'][0], vp2['bodies'][0]], ['PRS', f'{RRS_version}'], loc= 'best')
    plt.title(f'Distribution of features in PRS and {RRS_version}_1,_2,_3', fontsize= title_fontsize)
    # plt.grid(alpha= 0.2)
    plt.savefig(plot_path + f'/most_features_vp_PRS_{RRS_version}_1_2_3.pdf', bbox_inches= 'tight')
    print(f'Most features vp of {RRS_version} saved in {plot_path}.')
    plt.close()

    plt.figure(figsize= (6,6))
    plt.violinplot([-np.log10(PRS['Probability']), -np.log10(RRS['Probability'])], showmedians= True, widths= 0.8)
    plt.xticks([1,2], ['PRS', f'{RRS_version}'], fontsize= fontsize)
    plt.title(f'SLiM probabibility in PRS and {RRS_version}_1,_2,_3', fontsize= title_fontsize)
    plt.ylabel('SLiM Probability in -log10', fontsize= fontsize)
    plt.ylim([0, 16])
    # plt.grid(alpha= 0.2)
    plt.savefig(plot_path + f'/slim_probability_vp_PRS_{RRS_version}_1_2_3.pdf', bbox_inches= 'tight')
    print(f'SLiM probability vp of {RRS_version} saved in {plot_path}.')
    plt.close()

    plt.figure(figsize= (6,6))
    plt.violinplot([PRS['DomainFreqbyProtein'], RRS['DomainFreqbyProtein']], showmedians= True, widths= 0.8)
    plt.xticks([1,2], ['PRS', f'{RRS_version}'], fontsize= fontsize)
    plt.title(f'Domain frequency counted by protein in PRS and {RRS_version}_1,_2,_3', fontsize= title_fontsize)
    plt.ylabel('Domain frequency counted by protein', fontsize= fontsize)
    # plt.grid(alpha= 0.2)
    plt.savefig(plot_path + f'/domainfreqbyprotein_vp_PRS_{RRS_version}_1_2_3.pdf', bbox_inches= 'tight')
    print(f'DomainFreqbyProtein vp of {RRS_version} saved in {plot_path}.')
    plt.close()

    plt.figure(figsize= (6,6))
    plt.violinplot([PRS['DomainFreqinProteome'], RRS['DomainFreqinProteome']], showmedians= True, widths= 0.8)
    plt.xticks([1,2], ['PRS', f'{RRS_version}'], fontsize= fontsize)
    plt.title(f'Domain frequency in human proteome in PRS and {RRS_version}_1,_2,_3', fontsize= title_fontsize)
    plt.ylabel('Domain frequency in human proteome', fontsize= fontsize)
    # plt.grid(alpha= 0.2)
    plt.savefig(plot_path + f'/domainfreqinproteome_vp_PRS_{RRS_version}_1_2_3.pdf', bbox_inches= 'tight')
    print(f'DomainFreqinProteome vp of {RRS_version} saved in {plot_path}.')
    plt.close()

    # PRS has extreme values in DomainEnrichment_zscore. Make broken axis to show the distribution better.
    # N= 4 # portion for ax1
    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize= (6,6), sharey=True, gridspec_kw= {'wspace': 0.04, 'width_ratios': [N, 1]}) #gridspec adjusts space between two subplots
    #
    # # plot the same data on both axes
    # ax1.violinplot([PRS['DomainEnrichment_zscore'], RRS['DomainEnrichment_zscore']], vert= False)
    # ax2.violinplot([PRS['DomainEnrichment_zscore'], RRS['DomainEnrichment_zscore']], vert= False)
    #
    # # zoom-in / limit the view to different portions of the data
    # ax1.set_xlim(-5, 25)  # most of the data
    # ax2.set_xlim(85, 100)  # outliers only
    # ax1.grid(alpha= 0.2)
    # ax2.grid(alpha= 0.2)
    #
    # # hide the spines between ax and ax2
    # ax1.spines['right'].set_visible(False)
    # ax2.spines['left'].set_visible(False)
    # ax1.yaxis.tick_left()
    # ax1.tick_params(labelright=False)  # don't put tick labels at the top
    # ax2.yaxis.tick_right()
    #
    # # Now, let's turn towards the cut-out slanted lines.
    # # We create line objects in axes coordinates, in which (0,0), (0,1),
    # # (1,0), and (1,1) are the four corners of the axes.
    # # The slanted lines themselves are markers at those locations, such that the
    # # lines keep their angle and position, independent of the axes size or scale
    # # Finally, we need to disable clipping.
    #
    # d = .015  # how big to make the diagonal lines in axes coordinates
    # # arguments to pass to plot, just so we don't keep repeating them
    # kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    # ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)        # top-left diagonal
    # ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
    #
    # kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    # ax2.plot((-d * N, +d * N), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    # ax2.plot((-d * N, +d * N), (-d, +d), **kwargs)
    # plt.setp(ax2, yticks= [1,2], yticklabels= ['PRS','RRS'])
    #
    # fig.text(0.5, 0.04, 'Domain Enrichmnet z-score', ha='center', va='center')
    # fig.text(0.5, 0.92, f'Domain Enrichment z-score in PRS and {RRS_version}_1,_2,_3', ha='center', va='center')
    #
    # plt.savefig(plot_path + f'/domain_enrichment_zscore_vp_PRS_{RRS_version}_1_2_3.pdf', bbox_inches= 'tight')
    # print(f'Domain Enrichment z-score vp of {RRS_version} saved in {plot_path}.')
    # plt.close()
    N= 4
    f, (ax_top, ax_bottom) = plt.subplots(ncols=1, nrows=2, figsize= (6,6), sharex=True, gridspec_kw={'hspace':0.05, 'height_ratios': [1,N]})
    ax_top.violinplot([PRS['DomainEnrichment_zscore'], RRS['DomainEnrichment_zscore']])
    ax_bottom.violinplot([PRS['DomainEnrichment_zscore'], RRS['DomainEnrichment_zscore']])
    ax_top.set_ylim(85, 100)   # those limits are fake
    ax_bottom.set_ylim(-5, 25)

    ax_bottom.spines['top'].set_visible(False)
    ax_top.spines['bottom'].set_visible(False)
    ax_bottom.xaxis.tick_bottom()
    ax_top.xaxis.tick_top()

    ax = ax_top
    d = .015  # how big to make the diagonal lines in axes coordinates
    # arguments to pass to plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((-d, +d), (-d*N, +d*N), **kwargs)        # top-left diagonal

    ax2 = ax_bottom
    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal

    f.text(0.04, 0.5, 'Domain Enrichmnet z-score', ha='center', va='center', rotation= 'vertical', fontsize= fontsize)
    f.text(0.5, 0.92, f'Domain Enrichment z-score in PRS and {RRS_version}_1,_2,_3', ha='center', va='center', fontsize= fontsize)

    plt.setp(ax_bottom, xticks= [1.0, 2.0], xticklabels= ['PRS', 'RRS'])
    plt.xticks(fontsize= fontsize)
    plt.savefig(plot_path + f'/domain_enrichment_zscore_vp_PRS_{RRS_version}_1_2_3.pdf', bbox_inches= 'tight')
    print(f'Domain Enrichment z-score vp of {RRS_version} saved in {plot_path}.')
    plt.close()

def make_correlation_heatmap(PRS, RRS):
    """
    Plot the correlation matrix of all features and labels in the PRS and RRS. The plot is saved as pdf file.

    Args:
        PRS (pd.DataFrame): The processed PRS from the function preprocessing_dataset
        RRS (pd.DataFrame): The processed RRS (all triplicates concatenated into one RRS dataframe) from the function preprocessing_dataset
    """
    PRS['label']= 1.0
    RRS['label']= 0.0
    PRS_RRS_features= pd.concat([PRS, RRS], axis= 0)
    PRS_RRS_features= PRS_RRS_features[all_features_reordered_label]
    PRS_RRS_corr= PRS_RRS_features.corr()

    plt.figure(figsize= (8,8))
    ax= sns.heatmap(PRS_RRS_corr, vmin= -1, vmax= 1, center= 0, cmap= sns.diverging_palette(20, 220, n= 200), square= True, linewidths= 0.5, cbar_kws= {'shrink': 0.8})
    ax.set_xticklabels(ax.get_xticklabels(), rotation= 45, ha= 'right')
    ax.tick_params(axis= 'both', which= 'major', labelsize= fontsize - 1)
    ax.tick_params(axis= 'both', which= 'minor', labelsize= fontsize - 1)
    plt.title(f'Correlation heatmap of features and labels in PRS and {RRS_version}', fontsize= title_fontsize)
    plt.savefig(plot_path + f'/corr_heatmap_PRS_{RRS_version}_1_2_3.pdf', bbox_inches= 'tight')
    print(f'Correlation heatmap of {RRS_version} saved in {plot_path}.')
    plt.close()
    file_name= f'PRS_{RRS_version}_corr_table.tsv'
    PRS_RRS_corr.to_csv(plot_path[:-5] + file_name, sep= '\t')

if __name__ == '__main__':
    PRS = sys.argv[1]
    RRS_1= sys.argv[2]
    RRS_2= sys.argv[3]
    RRS_3= sys.argv[4]
    RRS_input_list= list([RRS_1, RRS_2, RRS_3])
    plot_path= '/'.join([i for i in RRS_1.split('/')[:-1]]) + '/Plots'
    RRS_version= RRS_1.split('/')[-2]

    PRS, RRS= preprocessing_dataset(PRS, RRS_input_list)
    make_DMI_fraction_plot(PRS, RRS)
    make_feature_violin_plots(PRS, RRS)
    make_correlation_heatmap(PRS, RRS)

    # python3 feature_analysis_plots_within_RRSv.py ../PRS/PRS_v3_only_human_with_pattern_alt_iso_swapped_removed_20210413_slim_domain_features_annotated.tsv ../RRS/RRSv1/RRSv1_1_20210427_slim_domain_features_annotated.tsv ../RRS/RRSv1/RRSv1_2_20210427_slim_domain_features_annotated.tsv ../RRS/RRSv1/RRSv1_3_20210427_slim_domain_features_annotated.tsv
    # python3 feature_analysis_plots_within_RRSv.py ../PRS/PRS_v3_only_human_with_pattern_alt_iso_swapped_removed_20210413_slim_domain_features_annotated.tsv ../RRS/RRSv2/RRSv2_1_20210428_slim_domain_features_annotated.tsv ../RRS/RRSv2/RRSv2_2_20210428_slim_domain_features_annotated.tsv ../RRS/RRSv2/RRSv2_3_20210428_slim_domain_features_annotated.tsv
    # python3 feature_analysis_plots_within_RRSv.py ../PRS/PRS_v3_only_human_with_pattern_alt_iso_swapped_removed_20210413_slim_domain_features_annotated.tsv ../RRS/RRSv3/RRSv3_1_20210428_slim_domain_features_annotated.tsv ../RRS/RRSv3/RRSv3_2_20210428_slim_domain_features_annotated.tsv ../RRS/RRSv3/RRSv3_3_20210428_slim_domain_features_annotated.tsv
    # python3 feature_analysis_plots_within_RRSv.py ../PRS/PRS_v3_only_human_with_pattern_alt_iso_swapped_removed_20210413_slim_domain_features_annotated.tsv ../RRS/RRSv4/RRSv4_1_20210428_slim_domain_features_annotated.tsv ../RRS/RRSv4/RRSv4_2_20210428_slim_domain_features_annotated.tsv ../RRS/RRSv4/RRSv4_3_20210428_slim_domain_features_annotated.tsv
