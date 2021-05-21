# Makes avg ROC and Precision Recall curve of all RRS versions in one plot
import sys
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

fontsize= 12
title_fontsize= 14
sns.set_style('darkgrid')
color= sns.color_palette('deep')

def make_stacked_PR_ROC_curves(RRS_paths): # input as list of RRS version paths

    # make ROC curve
    fig, ax= plt.subplots(figsize= (8,8))

    for c, RRS_path in zip(color, RRS_paths):
        RRS_version= RRS_path.split('/')[-2]
        with open(RRS_path + f'{RRS_version}_avg_ROC_scores.txt', 'r') as f:
            lines= [line.strip() for line in f.readlines()]
        mean_auc= float(lines[0].split(':')[1])
        mean_fpr= []
        mean_tpr= []
        std_tpr= []
        for line in lines[2:]:
            tabs= line.split('\t')
            mean_fpr.append(float(tabs[0]))
            mean_tpr.append(float(tabs[1]))
            std_tpr.append(float(tabs[2]))

        ax.plot(mean_fpr, mean_tpr, label= f'{RRS_version} mean AUC = {mean_auc:.2f}', color= c)
        ax.fill_between(mean_fpr, [(tpr - std) for tpr, std in zip(mean_tpr, std_tpr)], [(tpr + std) for tpr, std in zip(mean_tpr, std_tpr)], color= c, alpha= 0.3)
    ax.set(xlim=[-0.05, 1.05], ylim= [-0.05, 1.05])
    ax.legend(loc= 'lower right')
    ax.set_xlabel('False Positive Rate (Positive label: 1)', fontsize= fontsize)
    ax.set_ylabel('True Positive Rate (Positive label: 1)', fontsize= fontsize)
    ax.set_title(f'ROC curve averaged across triplicates of RRSv1,2,3,4 with ± 1 std', fontdict= {'fontsize': title_fontsize})

    # plt.grid(alpha= 0.2)
    plt.savefig(plot_path + f'/ROC_curve_RRSv1_2_3_4.pdf', bbox_inches= 'tight')
    print(f'ROC curve of RRSv1,2,3,4 saved in {plot_path}.')
    plt.close()

    # make PR curve
    fig, ax= plt.subplots(figsize= (8,8))

    for c, RRS_path in zip(color, RRS_paths):
        RRS_version= RRS_path.split('/')[-2]
        with open(RRS_path + f'{RRS_version}_avg_PR_scores.txt', 'r') as f:
            lines= [line.strip() for line in f.readlines()]
        mean_AP= float(lines[0].split(':')[1])
        mean_recall= []
        mean_precision= []
        std_precision= []
        for line in lines[2:]:
            tabs= line.split('\t')
            mean_recall.append(float(tabs[0]))
            mean_precision.append(float(tabs[1]))
            std_precision.append(float(tabs[2]))

        ax.plot(mean_recall, mean_precision, label= f'{RRS_version} mean AUC = {mean_auc:.2f}', color= c)
        ax.fill_between(mean_recall, [(prec - std) for prec, std in zip(mean_precision, std_precision)], [(prec + std) for prec, std in zip(mean_precision, std_precision)], color= c, alpha= 0.3)
    ax.set(xlim=[-0.05, 1.05], ylim= [0.45, 1.05])
    ax.legend(loc= 'lower left')
    ax.set_xlabel('Recall (Positive label: 1)', fontsize= fontsize)
    ax.set_ylabel('Precision (Positive label: 1)', fontsize= fontsize)
    ax.set_title(f'Precision Recall curve averaged across triplicates of RRSv1,2,3,4 with ± 1 std', fontdict= {'fontsize': title_fontsize})

    # plt.grid(alpha= 0.2)
    plt.savefig(plot_path + f'/PR_curve_RRSv1_2_3_4.pdf', bbox_inches= 'tight')
    print(f'PR curve of RRSv1,2,3,4 saved in {plot_path}.')
    plt.close()


if __name__ == '__main__':
    RRS_paths= list(sys.argv[1:])
    plot_path= '/'.join([i for i in sys.argv[1].split('/')[:-2]]) + '/Plots'

    make_stacked_PR_ROC_curves(RRS_paths)

# python3 make_stacked_PR_ROC_curves.py ../RRS/RRSv1/ ../RRS/RRSv2/ ../RRS/RRSv3/ ../RRS/RRSv4/
