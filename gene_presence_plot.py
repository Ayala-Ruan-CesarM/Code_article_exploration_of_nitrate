import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


df1=pd.read_csv("Nitrate_locus_per_taxa_filtered_manual.txt", sep="\t")
df2=pd.read_csv("Sulfate_locus_filtered.txt", sep="\t")

merged_df = df1.merge(df2, on="Taxa", how="outer") 
merged_df.to_csv("Test.txt", sep="\t", index=False)

data=pd.read_csv("Test.txt", sep="\t", index_col='Taxa')
df = data

# Create a list of columns to convert (all columns except the index)
columns_to_convert = df.columns
df[columns_to_convert] = df[columns_to_convert].replace(r'^\s*$', np.nan, regex=True)
df = df.astype(float, errors='ignore') # 'errors=ignore' handles non-numeric columns
df = df[df.sum(axis=1) != 1]
# A taxon is 'complete' if it has a non-zero, non-NA value for every gene.
list1=("narG", "narL", "narX", "narH", "narI", "narJ")
list2=("sat", "aprA", "aprB", "dsrA", "dsrB")

# Check for completeness based on list1
is_complete_list1 = (df[list(list1)] > 0).all(axis=1)

# Check for completeness based on list2
is_complete_list2 = (df[list(list2)] > 0).all(axis=1)

# A taxon is complete if it's complete in list1 OR complete in list2
complete_taxa = df[is_complete_list1 | is_complete_list2].index

#is_present = (df > 0) & (df.notna())
#complete_taxa = df[is_present.all(axis=1)].index

# --- 3. Reshape Data for Plotting (The "Tidy" Method) ---
# Melt the dataframe from wide to long format
plot_data = df.reset_index().rename(columns={'index': 'Taxa'}).melt(
    id_vars='Taxa',
    var_name='Gene',
    value_name='Value'
)

# Remove non-existent gene presences (zeros or NAs)
plot_data = plot_data.dropna(subset=['Value'])
plot_data = plot_data[plot_data['Value'] > 0]

# Add a 'Status' column to distinguish between complete and incomplete taxa
plot_data['Status'] = plot_data['Taxa'].apply(
    lambda x: 'Complete' if x in complete_taxa else 'Incomplete'
)

# --- 4. Create the Publication-Quality Dot Plot ---

# Set a professional theme and context for the plot
sns.set_theme(style="whitegrid", context="notebook")

# Create the figure and axes
fig, ax = plt.subplots(figsize=(12, 8))

# Define a color palette
palette = {"Complete": "#d62728", "Incomplete": "#1f77b4"} # Red and Blue

# Draw the scatterplot
sns.scatterplot(
    data=plot_data,
    x="Gene",
    y="Taxa",
    hue="Status",          # Color dots by the 'Status' column
    size="Value",          # Size dots by the 'Value' column
    sizes=(20, 200),       # Range of dot sizes
    palette=palette,       # Use our defined colors
    alpha=0.8,
    edgecolor="black",
    linewidth=1,
    ax=ax
)
ax.grid(False)
# --- 5. Customize and Refine the Plot ---

# Set title and labels with appropriate font sizes
ax.set_title('Gene Presence', fontsize=18, weight='bold', pad=20)
ax.set_xlabel('Genes', fontsize=14, weight='semibold')
ax.set_ylabel('Taxa', fontsize=14, weight='semibold')

# Customize tick labels
plt.xticks(rotation=45, ha='right', fontsize=10) # Rotate for readability
plt.yticks(fontsize=8, style='italic')

# Customize the legend
handles, labels = ax.get_legend_handles_labels()
# Separate legends for 'Status' (color) and 'Value' (size)
value_legend_index = labels.index('Value')
# Adjust the original legend to only show 'Value' (size)
value_legend = ax.legend(
    handles=handles[value_legend_index+1:],
    labels=[f'{float(l):.1f}' for l in labels[value_legend_index+1:]], # Format size labels
    title='Number of copies',
    bbox_to_anchor=(1.02, 0.7),
    loc='upper left',
    fontsize=10
)

# Final layout adjustments
plt.tight_layout(rect=[0, 0, 0.85, 1]) # Adjust layout to make space for legends
plt.savefig("plot3.pdf", dpi=1200)  # Perfect for publications or LaTeX
plt.show()
