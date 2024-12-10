'''
LOSO_revisit.py
'''
#
import os
import pandas as pd
import numpy as np
from scipy.stats.mstats import gmean

from sklearn.model_selection import LeaveOneGroupOut
from sklearn.linear_model import BayesianRidge, LinearRegression
from sklearn.model_selection import cross_val_score, RepeatedKFold
from sklearn.metrics import mean_absolute_error, r2_score
import seaborn as sns

from scipy.stats import spearmanr, pearsonr
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 16})
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# Functions

def clr(X):
#     return X ## DEBUG
    return np.log((X / gmean(X, axis=1).reshape(-1, 1)) )

# Calculate scores of the regression model
def calculateScore(true,pred):
    #mae, mape, r2, spearman, pearson

    overall_scores = {}
    overall_scores['mae'] = mean_absolute_error(true, pred)
    overall_scores['r2'] = r2_score(true, pred)
    overall_scores['spearman'] = spearmanr(true, pred)[0]
    overall_scores['pearson'] = pearsonr(true, pred)[0]
    # round to 3 digits
    for k,v in overall_scores.items():
        overall_scores[k] = round(v,3)
    return overall_scores


# turn scores dict to string for plot annotation
def scores2string(overall_scores):
    annotation_string = f"MAE = {overall_scores['mae']}\n" \
                        f"$R^2$ = {overall_scores['r2']}\n" \
                        f"Spearman = {overall_scores['spearman']}\n" \
                        f"Pearson = {overall_scores['pearson']}"
    return annotation_string


# Plot the regression
def plotRegResults(overall_scores,pred_dict, fig_size=(15,10), axis_title = 'train/validation', param_of_interest='Temperature',
                   fig_title = 'Bayesian Ridge - Temp prediction', DPI=300, margine_size=1,
                   do_save=False, saving_full_path=""):
    """
    :param scores_dict: {'mae','r2', 'spearman', 'pearson'}
    :param pred_dict: {'True','Pred','Group'}
    :return: figure
    """
    plt.rcParams.update({'font.size': 16})
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    # if nreads is a string that is a number, make it in a comma separated format (1,000 instead of 1000)

    annotation_string = scores2string(overall_scores)

    #y_val =np.array(pred_dict['True'])
    reported = np.array(pred_dict['True'])
    y_axis = reported.reshape(-1,1)
    #yval = y_val.reshape(-1,1)
    #preds = np.array(pred_dict['Pred'])
    inference = np.array(pred_dict['Pred'])
    x_axis = inference.reshape(-1,1)
    #yhat = preds.reshape(-1,1)
    # if TaraPolar is in a group, change it to Tara
    pred_dict['Group'] = ['Tara' if G=='TaraPolar' else G for G in pred_dict['Group']]

    match_color_dict = {'NEON':'#D28550','GEOTRACES':'#5FCAE6','Tara':'#123C67','Polarstern':'#9FD7C7', 'HOT-BATS':'#E78DC4'}
    targets = np.unique(pred_dict['Group'])
#     if len(targets)==1:
#         match_color_dict = {targets[0]:'#9FD7C7'}
#     else:
#         match_color_dict = {'NEON':'#D28550','GEOTRACES':'#5FCAE6','Tara':'#123C67'}
        
    # fig_title = f'Bayesian Ridge - Temp prediction\nk={k}, nreads={nreads}'
    fig = plt.figure(figsize = fig_size, dpi = DPI)
    fig.suptitle(fig_title,fontsize=18)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title(axis_title,fontsize=22)
    ax.set_xlabel(f'Inferred {param_of_interest}',fontsize=24)
    ax.set_ylabel(f'Reported {param_of_interest}',fontsize=24)
    
    for target, color in match_color_dict.items():
        indicesToKeep = [G==target for G in pred_dict['Group']]
        if sum(indicesToKeep)==0: continue
        #print(f'{target}: {sum(indicesToKeep)} points')
        ax.scatter(x = x_axis[indicesToKeep]
                   ,y =  y_axis[indicesToKeep]
                   , c = color
                   , marker = "D" if target=='NEON' else "o"
                   , s = 60, alpha=0.5
                   , label = fr'{target} ($n={sum(indicesToKeep)}$)')

    ax.tick_params(direction='inout')
    minimum_value = min(reported.min(),inference.min()) - margine_size
    maximum_value = max(reported.max(),inference.max()) + margine_size
    # identity line"
    ax.plot(range(int(minimum_value)-5,int(maximum_value)+5),range(int(minimum_value)-5,int(maximum_value)+5), c='k', alpha = 0.5, linestyle='--', label = 'Identity line')
    ax.set_xlim(minimum_value, maximum_value)
    ax.set_ylim(minimum_value, maximum_value)

    ax.plot(x_axis, LinearRegression().fit(x_axis, y_axis).predict(x_axis), c='k', label = 'Linear Regression')
    ax.grid(False)
    ax.legend(bbox_to_anchor=(1.04, 1),loc='upper left')

    ax.annotate(annotation_string,xy=(.02,.8),xycoords='axes fraction',
                fontsize=18, color = 'k',bbox=dict(boxstyle="round", fc="w"))
    #ax.legend(bbox_to_anchor=(1.05, 1), loc=2)
    if do_save:
        assert saving_full_path!="", "Provide full path to save the plot (saving_full_path) or change to do_save=False"
        plt.savefig(saving_full_path, format="pdf", bbox_inches="tight", transparent=True)
    return fig


# Parameters
is_rc = True
kmer = 4
is_mapped = True
param_of_interest = 'Temperature'

category = "Class" # Class or Phylum
tsv_directories = [
    f"/home/projects/zeevid/Analyses/2023-Tamir/NEON/Soil/all_samples_5k_freq/{category}",
    f"/home/projects/zeevid/Analyses/2023-Tamir/Geotraces/all_samples_5k_freq/{category}",
    f"/home/projects/zeevid/Analyses/2023-Tamir/Tara/PRJEB1787_prokaryote_DNA/all_samples_5k_freq/{category}"
]

# Initialize a dictionary to store the count of each Phylum/Class
valid_runs = pd.read_csv('/home/projects/zeevid/Analyses/2023-EbG/EnvironmentByGenome2024/mapping/valid_runs.csv')
NEON_naming_tuple = [(s,r.split('.fastq.gz;')[0].replace('_R1','').replace('-R1','')) for _,(s,r) in valid_runs.loc[valid_runs.project=='NEON',['sample','run_1']].iterrows() ]
valid_runs['RunID'] = ''
for s,r in NEON_naming_tuple:
    valid_runs.loc[valid_runs.loc[:,'sample']==s,'RunID'] = r

valid_runs.loc[valid_runs.project!='NEON','RunID'] = [x.split('_1')[0] for x in valid_runs.loc[valid_runs.project!='NEON','run_1']]

csv_df = pd.read_csv("/home/projects/zeevid/Analyses/2023-EbG/EnvironmentByGenome2024/data_recollection/Combined_datasets_full.csv")

work_df = valid_runs.copy().rename(columns={'sample':'UnifiedIndex'})
work_df = work_df.merge(csv_df[['UnifiedIndex','Temperature','LOSO_bySubsite']],on='UnifiedIndex',how='left')

category_counts = {}
for tsv_directory in tsv_directories:
    for file_name in os.listdir(tsv_directory):
        if file_name.endswith(".tsv") and file_name.split('.tsv')[0] in work_df.RunID.values:
            tsv_df = pd.read_csv(os.path.join(tsv_directory, file_name), sep='\t')
            elements_in_category = tsv_df.iloc[:, 0].tolist()
            for element in elements_in_category:
                if element in category_counts:
                    category_counts[element] += 1
                else:
                    category_counts[element] = 1

selected_elements_in_category = [element for element, count in category_counts.items() if count >= 100]

# Initialize lists to store later-on
mae_list = []
r2_list = []
pearson_list = []
spearman_list = []
dir_colors = ['dodgerblue', 'dodgerblue', 'mediumblue', 'forestgreen']
color_legend = ['Tara','TaraPolar', 'Geotraces','NEON']

feature_names = None  # Placeholder for feature names

dfs = {}
# Loop over each element in selected_elements_in_category (like phylum in phyla or class in classes)
for element in selected_elements_in_category:
    element_dict = {}
    dfs[element] = pd.DataFrame
    kmer_freq = []
    actual_temp = []
    dir_colors_list = []
    groups = []  # This will store the group labels for LOSO

    for tsv_directory in tsv_directories:
        dir_color = dir_colors[tsv_directories.index(tsv_directory)]
        for file_name in os.listdir(tsv_directory):
            if not (file_name.endswith(".tsv") and file_name.split('.tsv')[0] in work_df.RunID.values): continue
            temperature = work_df.loc[work_df.RunID==file_name.split('.tsv')[0],'Temperature'].values[0]

            tsv_df = pd.read_csv(os.path.join(tsv_directory, file_name), sep='\t')
            # if all values in element column are unique then set index. else raise error
            if not(tsv_df.iloc[:, 0].is_unique): raise ValueError("element column is not unique")
            tsv_df = tsv_df.set_index(tsv_df.columns[0])
            
            if feature_names is None:
                feature_names = tsv_df.columns.tolist()
            
            #element_df = tsv_df[tsv_df.iloc[:, 0] == element]
            if element not in tsv_df.index: continue
            name_key = work_df.loc[work_df.RunID==file_name.split('.tsv')[0],'UnifiedIndex'].values[0]
            element_row = tsv_df.loc[element]

            tetramer_freq = element_row.tolist()
            group = work_df.loc[work_df.RunID==file_name.split('.tsv')[0],'LOSO_bySubsite'].values[0]
            project = work_df.loc[work_df.RunID==file_name.split('.tsv')[0],'project'].values[0]
            element_dict[name_key] = [tetramer_freq, temperature, project, group]

    for i,key in enumerate(element_dict.keys()): # make a data frame for each element. Each row is a sample (name_key)
        
        temp_kmer_df = pd.DataFrame([element_dict[key][0]], columns=feature_names)
        temp_kmer_df['Temperature'] = element_dict[key][1]
        temp_kmer_df['Project'] = element_dict[key][2]
        temp_kmer_df['Group'] = element_dict[key][3]
        # row name is key
        temp_kmer_df.index = [key]
        if i==0:
            dfs[element] = temp_kmer_df.copy()
        else:
            dfs[element] = pd.concat([dfs[element],temp_kmer_df]).copy()
    
    remove_samples = dfs[element].index[(dfs[element]==0).any(axis=1)]
    dfs[element] = dfs[element].drop(remove_samples)

do_clr = True
cv_column = 'Group'

for element in selected_elements_in_category:
    #phylum = 'Actinobacteria'
    MainDf = dfs[element].copy()
    X_features = clr(MainDf[feature_names]) if do_clr else MainDf[feature_names]

    LOGOcv = LeaveOneGroupOut().split(\
                                    X = X_features.to_numpy(), \
                                    y = MainDf[param_of_interest].values, \
                                    groups = MainDf[cv_column].values \
                                    )

    n_iterations = MainDf[cv_column].nunique()
    CV_prediction_df = pd.DataFrame(columns=['ID','Actual','Pred','project'])
    all_scores = {}
    counter = 0
    validation_sizes=[]
    training_sizes=[]
    # cross-validation loop:
    for train_index, validation_index in LOGOcv:
        X_train, X_validation = X_features.to_numpy()[train_index],\
                                X_features.to_numpy()[validation_index]
        
        y_train, y_validation = MainDf[param_of_interest].values[train_index],\
                                MainDf[param_of_interest].values[validation_index]
        
        # we also collect the projects, for color coding
        projects_train, projects_validation = MainDf.Project.values[train_index],\
                                            MainDf.Project.values[validation_index]
        

        # choose your model. We work with bayesian ridge regression
        MLmodel = BayesianRidge().fit(X_train, y_train)
        scores = []
        scores.append(('train (R2)', MLmodel.score(X_train, y_train)))
        scores.append(('validation (MAE)', float(np.mean(np.abs(MLmodel.predict(X_validation)-y_validation) )) ))
        temp_pred = pd.DataFrame({'ID':MainDf.index.values[validation_index],
                                'Actual':y_validation,
                                'Pred':MLmodel.predict(X_validation),
                                'project':projects_validation})
        
        CV_prediction_df = pd.concat([CV_prediction_df,temp_pred])
        validation_sizes.append(len(validation_index))
        training_sizes.append(len(train_index))
        if counter%100==0:
            print(f"CV round {counter+1}/{n_iterations}:\n"
                f"#n_train: {len(train_index)}, #n_validation: {len(validation_index)}\n "
                f"{scores}\n")

        all_scores[counter] = scores
        counter += 1
    
    plotDF = CV_prediction_df.copy()

    ax_title = f"{param_of_interest} inference by {category} [{'CLR' if do_clr else 'NoCLR'}]"
    pred_dict = {'True':plotDF.Actual.values,
                'Pred':plotDF.Pred.values,
                'Group':plotDF.project.values,
                'ID':plotDF.ID.values}
    overall_scores = calculateScore(pred_dict['True'],pred_dict['Pred'])

    fig_cv = plotRegResults(overall_scores,pred_dict, fig_size=(10,10), axis_title = ax_title,param_of_interest=param_of_interest,
                            fig_title = f'LOSO CV - {element}', DPI=300, margine_size=2,
                                do_save=True, 
                                saving_full_path=f"/home/projects/zeevid/Analyses/2023-EbG/EnvironmentByGenome2024/GIT_IGNORED/Results/Taxonomy/Ridge_LOSO_{kmer}mers_{'rc' if is_rc else 'full'}_{category}_{element}.pdf")