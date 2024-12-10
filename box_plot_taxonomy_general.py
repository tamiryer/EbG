
# #works!! 
# import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
# import os

# # Load the data
# df = pd.read_csv('/home/projects/zeevid/tamirye/Data/Phyla_regression_with_outliers.csv')

# # Metrics to plot
# metrics = ['MAE_phyla', 'r2_phyla', 'Pearson_phyla', 'Spearman_phyla']
# bootstrap_metrics = ['bootstrap_MAE', 'bootstrap_r2', 'bootstrap_Pearson', 'bootstrap_Spearman']

# # Ensure output directory exists
# output_dir = '/home/projects/zeevid/tamirye/'
# os.makedirs(output_dir, exist_ok=True)

# # Loop over each metric to create separate figures
# for i, metric in enumerate(metrics):
#     plt.figure(figsize=(24, 18))
    
#     # Plot the boxplot for bootstrap data
#     ax = sns.boxplot(data=df[bootstrap_metrics[i]].apply(eval).tolist(), color='lightblue', width=0.5)

#        # Set transparency for the boxplot by adjusting the facecolor alpha
#     for patch in ax.artists:
#         patch.set_facecolor('lightblue')
#         patch.set_alpha(0.5)  # Set the transparency (0 is fully transparent, 1 is fully opaque)
        
#     # Adjust the thickness of the borders and the median line
#     for line in ax.lines:
#         line.set_linewidth(0.5)  # Set this to a smaller value to make lines thinner
    
    
#     # Overlay red dots for the phyla's singular value
#     sns.scatterplot(x=range(len(df)), y=df[metric], color='red', s=24, zorder=2, edgecolor='k', ax=ax)
    
#     # Set titles and labels
#     plt.title(metric.replace('_', ' '), fontsize=18) #.title()
#     plt.ylabel('Values', fontsize=18)
#     plt.xlabel('Phyla', fontsize=18)
#     plt.xticks(ticks=range(len(df)), labels=df['phylum'], rotation=45, fontsize=15)  # Adjust if there's a 'phylum_name' column
#     plt.yticks(fontsize=15)
#     # Adjust y-axis limits dynamically for better visibility
#     all_values = df[bootstrap_metrics[i]].apply(eval).explode().astype(float).values
#     min_value, max_value = min(all_values), max(all_values)
#     plt.ylim(min_value - 0.1, max_value + 0.1)  # Add some padding for better readability
    
#     plt.grid(True, linestyle='--', alpha=0.6)
#     plt.tight_layout()
    
#     # Save the plot as PDF
#     output_file = os.path.join(output_dir, f'boxplot_{metric}.pdf')
#     plt.savefig(output_file, format='pdf', dpi=300)
    
#     # Show the plot
#     plt.show()
    
#     print(f"Plot saved as {output_file}")


# import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
# import os

# # Load the data
# df = pd.read_csv('/home/projects/zeevid/tamirye/Data/Phyla_regression_with_outliers.csv')

# # Metrics to plot
# metrics = ['MAE_phyla', 'r2_phyla', 'Pearson_phyla', 'Spearman_phyla']
# bootstrap_metrics = ['bootstrap_MAE', 'bootstrap_r2', 'bootstrap_Pearson', 'bootstrap_Spearman']

# # Ensure output directory exists
# output_dir = '/home/projects/zeevid/tamirye/'
# os.makedirs(output_dir, exist_ok=True)

# # Loop over each metric to create separate figures
# for i, metric in enumerate(metrics):
#     plt.figure(figsize=(24, 16))  # Increase the figure size
    
#     # Plot the boxplot for bootstrap data
#     ax = sns.boxplot(data=df[bootstrap_metrics[i]].apply(eval).tolist(), color='lightblue', width=0.7)  # Increase the box width

#     # Set transparency for the boxplot by adjusting the facecolor alpha
#     for patch in ax.artists:
#         patch.set_facecolor('lightblue')
#         patch.set_alpha(0.5)  # Set the transparency (0 is fully transparent, 1 is fully opaque)
        
#     # Adjust the thickness of the borders and the median line
#     for line in ax.lines:
#         line.set_linewidth(1)  # Increase the line width

#     # Overlay red dots for the phyla's singular value
#     sns.scatterplot(x=range(len(df)), y=df[metric], color='red', s=50, zorder=2, edgecolor='k', linewidth=0.5)  # Increase the dot size and decrease the border width

#     # Set titles and labels
#     #plt.title(metric.replace('_', ' '), fontsize=20)
#     metric = metric.replace('_phyla', '')
#     plt.title(f'bootstrapping: {metric} comparison - phyla vs. general samples', fontsize=20)
#     plt.ylabel('Values', fontsize=20)
#     plt.xlabel('Phyla', fontsize=20)
#     plt.xticks(ticks=range(len(df)), labels=df['phylum'], rotation=45, fontsize=16)  # Adjust if there's a 'phylum_name' column
#     plt.yticks(fontsize=16)
    
#     # Adjust y-axis limits dynamically for better visibility
#     all_values = df[bootstrap_metrics[i]].apply(eval).explode().astype(float).values
#     min_value, max_value = min(all_values), max(all_values)
#     plt.ylim(min_value - 0.1, max_value + 0.1)  # Add some padding for better readability
    
#     plt.grid(True, linestyle='--', alpha=0.6)
#     plt.tight_layout()
    
#     # Save the plot as PDF
#     output_file = os.path.join(output_dir, f'boxplot_{metric}_phyla.pdf')
#     plt.savefig(output_file, format='pdf', dpi=300)
    
#     # Show the plot
#     plt.show()
    
#     print(f"Plot saved as {output_file}")

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# Load the data
df = pd.read_csv('/home/projects/zeevid/tamirye/Data/Phyla_regression_with_outliers.csv')

# Metrics to plot
metric = 'Spearman_phyla'
bootstrap_metric = 'bootstrap_Spearman'

# Ensure output directory exists
output_dir = '/home/projects/zeevid/tamirye/'
os.makedirs(output_dir, exist_ok=True)

plt.figure(figsize=(18, 12))

# Plot the boxplot for bootstrap data
ax = sns.boxplot(data=df[bootstrap_metric].apply(eval).tolist(), color='lightblue', width=0.7)

# Set transparency for the boxplot by adjusting the facecolor alpha
for patch in ax.artists:
    patch.set_facecolor('lightblue')
    patch.set_alpha(0.5)

# Adjust the thickness of the borders and the median line
for line in ax.lines:
    line.set_linewidth(1)

# Overlay red dots for the phyla's singular value
sns.scatterplot(x=range(len(df)), y=df[metric], color='red', s=50, zorder=2, edgecolor='k', linewidth=0.5)

# Set titles and labels
metric = metric.replace('_phyla', '')
plt.title(f'Bootstrapping: {metric} comparison - Phyla vs. general samples', fontsize=20)
plt.ylabel('Values', fontsize=20)
plt.xlabel('Phyla', fontsize=20)
plt.xticks(ticks=range(len(df)), labels=df['phylum'], rotation=60, fontsize=17, ha='right') # Align text to the right)
plt.yticks(fontsize=16)

# Adjust y-axis limits to start from 0.7
plt.ylim(0.70, 1.0)

plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()

# Save the plot as PDF
output_file = os.path.join(output_dir, f'boxplot_{metric}_phyla.png')
plt.savefig(output_file, format='png', dpi=300)

# Show the plot
plt.show()

print(f"Plot saved as {output_file}")