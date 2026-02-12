import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches 


# Data from the provided text, manually adjusted to match the boxplot's displayed values
# and to include a 'Treatment Name' for better labeling as seen in the boxplot.
# Values greater than 22.5 (the visible limit in the example plot) are capped at 22.5
# for visualization purposes in the original plot, but here we'll use the original values
# and extend the y-axis.

data = {
    'Tratamiento': [
        'Control', 'Control', 'Control',
        'Nitrate_HC', 'Nitrate_HC', 'Nitrate_HC',
        'Nitrate_MC', 'Nitrate_MC', 'Nitrate_MC',
        'Nitrate_LC', 'Nitrate_LC', 'Nitrate_LC',
        'Nitrite_HC','Nitrite_HC','Nitrite_HC',
        'Nitrite_MC','Nitrite_MC','Nitrite_MC',
        'Nitrite_LC','Nitrite_LC','Nitrite_LC',
        'Molibdate_HC', 'Molibdate_HC', 'Molibdate_HC',
        'Molibdate_MC', 'Molibdate_MC', 'Molibdate_MC',        
        'Molibdate_LC', 'Molibdate_LC', 'Molibdate_LC',
        'Selenite _HC', 'Selenite _HC', 'Selenite _HC',
        'Selenite _MC', 'Selenite _MC', 'Selenite _MC',
        'Selenite _LC', 'Selenite _LC', 'Selenite _LC',
        'Nit_Molib_HC', 'Nit_Molib_HC', 'Nit_Molib_HC',
        'Nit_Molib_MC', 'Nit_Molib_MC', 'Nit_Molib_MC',
        'Nit_Molib_LC', 'Nit_Molib_LC', 'Nit_Molib_LC',
        'Nit_Sele_HC', 'Nit_Sele_HC', 'Nit_Sele_HC',
        'Nit_Sele_MC', 'Nit_Sele_MC', 'Nit_Sele_MC',
        'Nit_Sele_LC', 'Nit_Sele_LC', 'Nit_Sele_LC',
        'Nit_Nitr_HC', 'Nit_Nitr_HC', 'Nit_Nitr_HC',
        'Nit_Nitr_MC', 'Nit_Nitr_MC', 'Nit_Nitr_MC',
        'Nit_Nitr_LC', 'Nit_Nitr_LC', 'Nit_Nitr_LC',
    ],
    'Numero de Días': [
        4, 5, 5,  # Control 
        28, 28, 28,  # Tratamiento 1 - Alta (Nitrate_HC) 
        28, 28, 28,  # Tratamiento 1 - Media (Nitrate_MC) - This is guessed to be 'MC' from boxplot, assuming "Media" is Medium
        28, 28, 28, # Tratamiento 1 - Baja (Nitrate_LC) - This is guessed to be 'LC' from boxplot, assuming "Baja" is Low
        28, 28, 28,  # Tratamiento 2 - Alta (Nitrate_HC) 
        28, 28, 28,  # Tratamiento 2 - Media (Nitrate_MC) - This is guessed to be 'MC' from boxplot, assuming "Media" is Medium
        28, 28, 28,  # Tratamiento 2
        11, 11, 11,  # Tratamiento 3 - Alta (Molibdate_HC) 
        7, 7, 7,  # Tratamiento 3 - Media (Molibdate_MC) 
        4, 4, 4,  # Tratamiento 3 - Baja (Molibdate_LC) 
        15, 19, 15,  # Tratamiento 4 - Alta (Selenite _HC) 
        7, 5, 7,  # Tratamiento 4 - Media (Selenite _MC) 
        4, 4, 4,  # Tratamiento 4 - Baja (Selenite _LC) 
        28, 28, 28,  # Tratamiento 5 - Alta (Nt_Molib_HC) 
        10, 28, 28,  # Tratamiento 5 - Media (Nt_Molib_MC) 
        28, 28, 28,  # Tratamiento 5 - Baja (Nt_Molib_LC) 
        22, 22, 4,  # Tratamiento 6 - Alta (Nt_Sele_HC) 
        13, 15, 22,  # Tratamiento 6 - Media (Nt_Sele_MC) 
        11, 7, 10,  # Tratamiento 6 - Baja (Nt_Sele_LC) 
        28, 28, 28,  # Tratamiento 7 - Alta (Nt_Nitr_HC) 
        28, 28, 28,  # Tratamiento 7 - Media (Nt_Nitr_MC) 
        28, 28, 28,  # Tratamiento 7 - Baja (Nt_Nitr_LC)  
    ]
}

df = pd.DataFrame(data)

treatment_colors = {
    'Control': 'red',
    'Nitrate_HC': 'skyblue',
    'Nitrate_MC': 'lightgreen',
    'Nitrate_LC': 'lightcoral',
    'Nitrite_HC': 'chocolate',
    'Nitrite_MC': 'violet',
    'Nitrite_LC': 'darkred',
    'Molibdate_HC': 'gold',
    'Molibdate_MC': 'purple',
    'Molibdate_LC': 'orange',
    'Selenite _HC': 'darkcyan',
    'Selenite _MC': 'brown',    
    'Selenite _LC': 'magenta',
    'Nit_Molib_HC': 'darkblue',
    'Nit_Molib_MC': 'gray',
    'Nit_Molib_LC': 'pink',
    'Nit_Sele_HC': 'teal',
    'Nit_Sele_MC': 'maroon',
    'Nit_Sele_LC': 'lime',
    'Nit_Nitr_HC': 'indigo',
    'Nit_Nitr_MC': 'darkgreen',
    'Nit_Nitr_LC': 'olive',
}

# Set the style of the plots
sns.set_style("whitegrid")

plt.figure(figsize=(15, 8)) # Adjust figure size as needed

# Create the swarm plot with specific colors and a legend
# Assign the Axes object to a variable
ax = sns.swarmplot(
    x='Tratamiento',
    y='Numero de Días',
    data=df,
    hue='Tratamiento', # Use 'Tratamiento' to differentiate colors and create legend entries
    palette=treatment_colors, # Apply the custom color palette
    size=10,
    dodge=False # This helps separate points from different treatments if they overlap on the x-axis
)

# Improve the Y-axis: Extend beyond 28 days
ax.set_ylim(0, 30) # Set the y-axis limit from 0 to 30 to clearly see values up to 28 and beyond

# --- NEW: Adjust X-axis limits for padding ---
# Get the number of unique treatments to calculate appropriate x-axis limits
num_treatments = len(df['Tratamiento'].unique())
# Set x-axis limits with a small padding (e.g., 0.5 units on each side)
ax.set_xlim(-0.9, num_treatments - 0.1)
# --- NEW: Remove X-axis labels and ticks ---
# This line removes the numerical tick marks that matplotlib automatically places.
ax.set_xticks([])
# This line removes the text labels that correspond to the tick marks (your treatment names).
ax.set_xticklabels([])
# This line removes the overall label for the x-axis, as it's no longer necessary without individual treatment labels.
ax.set_xlabel('')
ax.set_title('Treatment Comparison', fontsize=16)
#ax.set_xlabel('Tratamientos', fontsize=12)
ax.set_ylabel('Days of incubation', fontsize=12)
plt.xticks(rotation=90) # Rotate x-axis labels for better readability

# Manually create legend handles and labels
legend_handles = [mpatches.Patch(color=treatment_colors[label], label=label)
                  for label in treatment_colors.keys()]

# Create the legend explicitly using the collected handles and labels
ax.legend(handles=legend_handles, title='Treatment', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

plt.tight_layout(rect=[0, 0, 0.90, 1]) # Adjust layout to make space for the legend

plt.show()
