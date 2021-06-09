import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, plot_confusion_matrix, plot_roc_curve, plot_precision_recall_curve, f1_score
from sklearn.model_selection import cross_val_score
from sklearn.tree import plot_tree
from joblib import dump, load

all_features= ['Probability', 'IUPredShort', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLC', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'DomainFreqbyProtein1', 'DomainFreqinProteome1']
all_features_renamed= ['Probability', 'IUPredShort', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLC', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'DomainFreqbyProtein', 'DomainFreqinProteome']

fontsize= 12
title_fontsize= 14
sns.set_style('darkgrid')

def preprocessing_dataset(PRS_input, RRS_input): # takes the PRS and RRS, concatenate them and preprocessing the NaNs and dummy value.
    PRS= pd.read_csv(PRS_input, sep= '\t', index_col= 0)
    PRS['label']= 1
    RRS= pd.read_csv(RRS_input, sep= '\t', index_col= 0)
    RRS['label']= 0
    for df in [PRS, RRS]:
        df.replace(88888, df.DomainEnrichment_zscore.median(), inplace= True)
        for ind, row in df.iterrows():
            if pd.notna(row['DomainFreqbyProtein2']):
                df.loc[ind, 'DomainFreqbyProtein1'] = np.mean([row['DomainFreqbyProtein1'], row['DomainFreqbyProtein2']])
                df.loc[ind, 'DomainFreqinProteome1'] = np.mean([row['DomainFreqinProteome1'], row['DomainFreqinProteome2']])
    df= pd.concat([PRS, RRS], axis= 0, ignore_index= True)
    df.dropna(subset= all_features, inplace= True)
    df.rename(columns= {'DomainFreqbyProtein1': 'DomainFreqbyProtein', 'DomainFreqinProteome1': 'DomainFreqinProteome'}, inplace= True)
    X= df[all_features_renamed].copy()
    y= df['label']
    return df, X, y

def split_fit_rf(X, y, exclude_feature= None): # features to be dropped saved as list
    rf= RandomForestClassifier(n_estimators= 1000, oob_score= True, verbose= True, n_jobs= -1)

    if exclude_feature != None:
        X= X.drop(labels= exclude_feature, axis= 1)
    X_train, X_test, y_train, y_test= train_test_split(X, y, random_state= 0, stratify= y)
    rf.fit(X_train, y_train)

    return rf, X_train, X_test, y_train, y_test

def make_confusion_matrix(split_fit_outputs):
    fig, axes= plt.subplots(1,3, figsize= (12, 12))

    for i, ele in enumerate(zip(split_fit_outputs, axes)):
        output, ax= ele
        rf, X_train, X_test, y_train, y_test= output
        plot_confusion_matrix(rf, X_test, y_test, ax= ax, colorbar= False)
        ax.set_title(f'RRSv1_{i+1}', fontsize= fontsize)

    plt.savefig(plot_path + f'/confusion_matrix_{RRS_version}_1_2_3.pdf', bbox_inches= 'tight')
    print(f'Confusion matrix of {RRS_version} saved in {plot_path}.')
    plt.close()

def make_ROC_curve(split_fit_outputs):

    tprs= []
    aucs= []
    mean_fpr= np.linspace(0, 1, 100)

    fig, ax= plt.subplots(figsize= (8,8))

    for i, output in enumerate(split_fit_outputs):
        rf, X_train, X_test, y_train, y_test= output
        roc= plot_roc_curve(rf, X_test, y_test, ax= ax, name= f'{RRS_version}_{i+1}')
        interp_tpr= np.interp(mean_fpr, roc.fpr, roc.tpr)
        interp_tpr[0]= 0.0
        tprs.append(interp_tpr)
        aucs.append(roc.roc_auc)

    # plt.grid(alpha= 0.2)
    plt.title(f'Receiver Operating Characteristics curve ({RRS_version})')
    plt.savefig(plot_path + f'/ROC_curve_{RRS_version}_1_2_3.pdf', bbox_inches= 'tight')
    print(f'ROC curve of {RRS_version} saved in {plot_path}.')
    plt.close()

    # plot avg ROC curve across triplicates
    mean_tpr= np.mean(tprs, axis= 0)
    mean_tpr[-1]= 1.0
    std_tpr= np.std(tprs, axis= 0)
    mean_auc= np.mean(aucs)

    fig, ax= plt.subplots(figsize= (8,8))

    ax.plot(mean_fpr, mean_tpr, label= f'{RRS_version} mean AUC = {mean_auc:.2f}')
    ax.fill_between(mean_fpr, mean_tpr - std_tpr, mean_tpr + std_tpr, alpha= 0.3, label= '1 std. dev.')
    ax.set(xlim=[-0.05, 1.05], ylim= [-0.05, 1.05])
    ax.legend(loc= 'lower right')
    ax.set_xlabel('False Positive Rate (Positive label: 1)', fontsize= fontsize)
    ax.set_ylabel('True Positive Rate (Positive label: 1)', fontsize= fontsize)
    ax.set_title(f'Receiver Operating Characteristics curve averaged across triplicates of {RRS_version}', fontdict= {'fontsize': title_fontsize})

    # plt.grid(alpha= 0.2)
    plt.savefig(plot_path + f'/avg_ROC_curve_{RRS_version}_1_2_3.pdf', bbox_inches= 'tight')
    print(f'Average ROC curve of {RRS_version} saved in {plot_path}.')
    plt.close()

    # write avg ROC into a txt file so that it can be plotted with other RRS versions
    with open(plot_path[:-5] + f'{RRS_version}_avg_ROC_scores.txt', 'w') as f:
        f.write(f'mean_auc:{mean_auc:.2f}\n')
        f.write('mean_fpr\tmean_tpr\tstd_tpr\n')
        for ele in zip(mean_fpr, mean_tpr, std_tpr):
            f.write(f'{ele[0]:.4f}\t{ele[1]:.4f}\t{ele[2]:.4f}')
            f.write('\n')
    print(f'ROC scores saved in {plot_path[:-5]}{RRS_version}_avg_ROC_scores.txt.')

def make_precision_recall_curve(split_fit_outputs):

    precisions= []
    aps= []
    mean_recall= np.linspace(0.0, 1.0, 100)

    fig, ax= plt.subplots(figsize= (8,8))

    for i, output in enumerate(split_fit_outputs):
        rf, X_train, X_test, y_train, y_test= output
        prc= plot_precision_recall_curve(rf, X_test, y_test, ax= ax, name= f'{RRS_version}_{i+1}')
        interp_prec= np.interp(mean_recall, np.flipud(prc.recall), np.flipud(prc.precision))
        interp_prec[0]= 1.0
        interp_prec[-1]= len(y_test[y_test == 1])/ len(y_test)
        precisions.append(interp_prec)
        aps.append(prc.average_precision)

    # plt.grid(alpha= 0.2)
    plt.title(f'Precision Recall curve ({RRS_version})')
    plt.savefig(plot_path + f'/precision_recall_curve_{RRS_version}_1_2_3.pdf', bbox_inches= 'tight')
    print(f'Precision recall curve of {RRS_version} saved in {plot_path}.')
    plt.close()

    # plot avg PR curve across triplicates
    mean_precision= np.mean(precisions, axis= 0)
    # mean_precision[-1]= 0.5
    std_precision= np.std(precisions, axis= 0)
    mean_ap= np.mean(aps)

    fig, ax= plt.subplots(figsize= (8,8))

    ax.plot(mean_recall, mean_precision, label= f'{RRS_version} mean AP = {mean_ap:.2f}')
    ax.fill_between(mean_recall, mean_precision - std_precision, mean_precision + std_precision, alpha= 0.3, label= '1 std. dev.')
    ax.set(xlim=[-0.05, 1.05], ylim= [0.45, 1.05])
    ax.legend(loc= 'lower left')
    ax.set_xlabel('Recall (Positive label: 1)', fontsize= fontsize)
    ax.set_ylabel('Precision (Positive label: 1)', fontsize= fontsize)
    ax.set_title(f'Precision Recall curve averaged across triplicates of {RRS_version}', fontdict= {'fontsize': title_fontsize})

    # plt.grid(alpha= 0.2)
    plt.savefig(plot_path + f'/avg_PR_curve_{RRS_version}_1_2_3.pdf', bbox_inches= 'tight')
    print(f'Average PR curve of {RRS_version} saved in {plot_path}.')
    plt.close()

    # write avg PR into a txt file so that it can be plotted with other RRS versions
    with open(plot_path[:-5] + f'{RRS_version}_avg_PR_scores.txt', 'w') as f:
        f.write(f'mean_ap:{mean_ap:.2f}\n')
        f.write('mean_recall\tmean_precision\tstd_precision\n')
        for ele in zip(mean_recall, mean_precision, std_precision):
            f.write(f'{ele[0]:.4f}\t{ele[1]:.4f}\t{ele[2]:.4f}')
            f.write('\n')
    print(f'PR scores saved in {plot_path[:-5]}{RRS_version}_avg_PR_scores.txt.')

def make_cvacc_oob_acc_plot(split_fit_outputs):

    cvacc= []
    cvacc_std= []
    oob= []
    acc= []

    for i, output in enumerate(split_fit_outputs):
        rf, X_train, X_test, y_train, y_test= output
        oob.append(rf.oob_score_)
        acc.append(accuracy_score(y_test, rf.predict(X_test)))
        cvacc.append(np.mean(cross_val_score(RandomForestClassifier(n_estimators= 1000), X_train, y_train, cv= 5, n_jobs= -1)))
        cvacc_std.append(np.std(cross_val_score(RandomForestClassifier(n_estimators= 1000), X_train, y_train, cv= 5, n_jobs= -1)))

    N= 3
    ind= np.arange(N)
    width= 0.25

    plt.figure(figsize= (8,6))

    plt.bar(ind, cvacc, width, yerr= cvacc_std, capsize= 7, color= 'c', label= '5-fold CV')
    plt.bar(ind + width, oob, width, color= 'b', label= 'oob')
    plt.bar(ind + 2 * width, acc, width, color= 'k', label= 'holdout')

    plt.xticks(ind + width, [f'{RRS_version}_1', f'{RRS_version}_2', f'{RRS_version}_3'], fontsize= fontsize)
    plt.ylim(0, 1)
    plt.title(f'Comparison of accuracy scores computed with different strategies ({RRS_version})', fontsize= title_fontsize)
    plt.ylabel('Accuracy score', fontsize= fontsize)
    plt.legend(bbox_to_anchor= (1.05, 1), loc= 'upper left', borderaxespad= 0.)
    # plt.grid(alpha= 0.2)

    plt.savefig(plot_path + f'/accuracy_scores_{RRS_version}_1_2_3.pdf', bbox_inches= 'tight')
    print(f'Acurracy scores of {RRS_version} saved in {plot_path}.')
    plt.close()

def make_feature_importance_plot(split_fit_outputs, exclude_feature= None):
    for i, output in enumerate(split_fit_outputs):
        rf, _, _, _, _= output
        if i== 0:
            if exclude_feature != None:
                features= list(filter(lambda feat: feat not in exclude_feature, all_features_renamed))
                xlim= None
            else:
                features= all_features_renamed
                xlim= [0, 0.16]
            feat_imp_df= pd.DataFrame(data= {'Features': features, f'{RRS_version}_{i+1}': rf.feature_importances_})
        else:
            feat_imp_df= pd.concat([feat_imp_df, pd.Series(rf.feature_importances_, name= f'{RRS_version}_{i+1}')], axis= 1)
    feat_imp_df['mean']= feat_imp_df.mean(axis= 1)
    feat_imp_df['std']= feat_imp_df.std(axis= 1)
    feat_imp_df= feat_imp_df.sort_values(by= f'{RRS_version}_1', ascending= True)

    N= len(features)
    ind= np.arange(N)
    height= 0.3

    plt.figure(figsize= (8,10))
    plt.barh(ind + 2*height, feat_imp_df[f'{RRS_version}_1'], height, color= 'c', label= f'{RRS_version}_1')
    plt.barh(ind + height, feat_imp_df[f'{RRS_version}_2'], height, color= 'b', label= f'{RRS_version}_2')
    plt.barh(ind, feat_imp_df[f'{RRS_version}_3'], height, color= 'r', label= f'{RRS_version}_3')

    plt.yticks(ind + height, feat_imp_df.Features)
    plt.legend(loc= 'best')
    plt.title(f'Feature importance of {RRS_version}')
    plt.xlabel('Gini index')
    # plt.grid(alpha= 0.2)

    plt.savefig(plot_path + f'/feature_importance_{RRS_version}_1_2_3.pdf', bbox_inches= 'tight')
    print(f'Feature importance plot of {RRS_version} saved in {plot_path}.')
    plt.close()

    # plot avg and std of feature importance across RRS triplicates.
    feat_imp_df= feat_imp_df.sort_values(by= 'mean', ascending= True)

    plt.figure(figsize= (6, 8))

    plt.barh(feat_imp_df['Features'], feat_imp_df['mean'], alpha= 0.9, color= 'steelblue', xerr= feat_imp_df['std'], capsize= 5)

    plt.xlim(xlim)
    plt.title(f'Average of feature importance across triplicates of {RRS_version}', fontsize= fontsize)
    plt.xlabel('Gini index', fontsize= fontsize)
    # plt.grid(alpha= 0.2)

    plt.savefig(plot_path + f'/avg_std_feature_importance_{RRS_version}_1_2_3.pdf', bbox_inches= 'tight')
    print(f'Average feature importance plot of {RRS_version} saved in {plot_path}.')
    plt.close()

def plot_one_tree(split_fit_outputs):
    for i, output in enumerate(split_fit_outputs):
        rf, _, _, _, _= output
        plt.figure(figsize= (16,16))
        plot_tree(rf.estimators_[0], filled= True, proportion= True, rounded= True, fontsize= 5, feature_names= all_features_renamed, precision= 2)
        plt.savefig(plot_path + f'/plot_first_tree_{RRS_version}_{i+1}.pdf', bbox_inches= 'tight')
        print(f'First tree of {RRS_version}_{i+1} saved in {plot_path}.')
        plt.close()

def write_classification_report(split_fit_outputs, exclude_feature= None):

    if exclude_feature != None:
        if len(exclude_feature) > 10:
            target_name= plot_path[:-5] + f'RF_classification_report_{RRS_version}_with_{"_".join(feat for feat in all_features_renamed if feat not in exclude_feature)}.txt'
        else:
            target_name= plot_path[:-5] + f'RF_classification_report_{RRS_version}_no_{"_".join(feat for feat in exclude_feature)}.txt'
    else:
        target_name= plot_path[:-5] + f'RF_classification_report_{RRS_version}.txt'
    with open(target_name, 'w') as f:
        for i, output in enumerate(split_fit_outputs):
            rf, X_train, X_test, y_train, y_test= output
            f.write(f'{RRS_version}_{i+1}\n')
            if exclude_feature != None:
                f.write(f'Trained on {rf.n_features_} features: {",".join([feat for feat in all_features_renamed if feat not in exclude_feature])}.\n')
            else:
                f.write(f'Trained on all {rf.n_features_} features.\n')
            f.write(classification_report(y_test, rf.predict(X_test)))
            f.write('\n')
    print(f'Classification report on {RRS_version} saved in {plot_path[:-5]}.')

def write_prediction(split_fit_outputs, df_list):

    for i, ele in enumerate(zip(split_fit_outputs, df_list)):
        output, df= ele
        rf, X_train, X_test, y_train, y_test= output
        test_index= y_test.index
        test_df= df.loc[test_index]
        test_df['predicted_values']= rf.predict(X_test)
        test_df['predicted_prob_1']= rf.predict_proba(X_test)[:, 1]
        # df= pd.DataFrame(data= {'actual_values': y_test, 'predicted_values': rf.predict(X_test), 'predicted_prob_1': rf.predict_proba(X_test)[:, 1]})
        test_df.to_csv(plot_path[:-5] + f'RF_scoring_on_test_set_{RRS_version}_{i+1}.tsv', sep= '\t')

def save_RF(split_fit_outputs):

    for i, output in enumerate(split_fit_outputs):
        rf, _, _, _, _ = output
        with open(plot_path[:-5] + f'RF_model_{RRS_version}_{i+1}.joblib', 'wb') as f:
            dump(rf, f)
            # load by executing rf= load(f)
        print(f'Fitted RandomForestClassifier with {RRS_version}_{i+1} saved in {plot_path[:-5]}.')

if __name__ == '__main__':
    PRS= sys.argv[1]
    RRS_1= sys.argv[2]
    RRS_2= sys.argv[3]
    RRS_3= sys.argv[4]
    plot_path= '/'.join([i for i in RRS_1.split('/')[:-1]]) + '/Plots'
    RRS_version= RRS_1.split('/')[-2]

    df_1, X_1, y_1= preprocessing_dataset(PRS, RRS_1)
    df_2, X_2, y_2= preprocessing_dataset(PRS, RRS_2)
    df_3, X_3, y_3= preprocessing_dataset(PRS, RRS_3)

    output_1= split_fit_rf(X_1, y_1)
    output_2= split_fit_rf(X_2, y_2)
    output_3= split_fit_rf(X_3, y_3)
    outputs= [output_1, output_2, output_3]
    make_confusion_matrix(outputs)
    make_ROC_curve(outputs)
    make_precision_recall_curve(outputs)
    make_cvacc_oob_acc_plot(outputs)
    make_feature_importance_plot(outputs)
    plot_one_tree(outputs)
    write_classification_report(outputs)
    write_prediction(outputs, [df_1, df_2, df_3])
    save_RF(outputs)

    # train model without the least important feature
    # exclude_feature= ['IUPredLong', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLC', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'DomainFreqbyProtein', 'DomainFreqinProteome']
    # output_1= split_fit_rf(X_1, y_1, exclude_feature= exclude_feature)
    # output_2= split_fit_rf(X_2, y_2, exclude_feature= exclude_feature)
    # output_3= split_fit_rf(X_3, y_3, exclude_feature= exclude_feature)
    # outputs= [output_1, output_2, output_3]
    # make_feature_importance_plot(outputs, exclude_feature= exclude_feature)
    # write_classification_report(outputs, exclude_feature= exclude_feature)

# python3 RandomForest_fitting.py ../PRS/PRS_v3_only_human_with_pattern_alt_iso_swapped_removed_20210413_slim_domain_features_annotated.tsv ../RRS/RRSv1/RRSv1_1_20210427_slim_domain_features_annotated.tsv ../RRS/RRSv1/RRSv1_2_20210427_slim_domain_features_annotated.tsv ../RRS/RRSv1/RRSv1_3_20210427_slim_domain_features_annotated.tsv
# python3 RandomForest_fitting.py ../PRS/PRS_v3_only_human_with_pattern_alt_iso_swapped_removed_20210413_slim_domain_features_annotated.tsv ../RRS/RRSv2/RRSv2_1_20210428_slim_domain_features_annotated.tsv ../RRS/RRSv2/RRSv2_2_20210428_slim_domain_features_annotated.tsv ../RRS/RRSv2/RRSv2_3_20210428_slim_domain_features_annotated.tsv
# python3 RandomForest_fitting.py ../PRS/PRS_v3_only_human_with_pattern_alt_iso_swapped_removed_20210413_slim_domain_features_annotated.tsv ../RRS/RRSv3/RRSv3_1_20210428_slim_domain_features_annotated.tsv ../RRS/RRSv3/RRSv3_2_20210428_slim_domain_features_annotated.tsv ../RRS/RRSv3/RRSv3_3_20210428_slim_domain_features_annotated.tsv
# python3 RandomForest_fitting.py ../PRS/PRS_v3_only_human_with_pattern_alt_iso_swapped_removed_20210413_slim_domain_features_annotated.tsv ../RRS/RRSv4/RRSv4_1_20210428_slim_domain_features_annotated.tsv ../RRS/RRSv4/RRSv4_2_20210428_slim_domain_features_annotated.tsv ../RRS/RRSv4/RRSv4_3_20210428_slim_domain_features_annotated.tsv
