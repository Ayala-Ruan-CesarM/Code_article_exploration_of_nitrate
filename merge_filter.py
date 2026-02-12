import pandas as pd

##############################Nitrate##################################
# Esto une los dataframes 
df1=pd.read_csv("narG_per_specie.txt", sep="\t")
df2=pd.read_csv("narL_per_specie.txt", sep="\t")
df3=pd.read_csv("narX_per_specie.txt", sep="\t")
df4=pd.read_csv("narH_per_specie.txt", sep="\t")
df5=pd.read_csv("narI_per_specie.txt", sep="\t")
df6=pd.read_csv("narJ_per_specie.txt", sep="\t")

merged_df = df1.merge(df2, on="Taxa", how="outer") \
               .merge(df3, on="Taxa", how="outer") \
               .merge(df4, on="Taxa", how="outer") \
               .merge(df5, on="Taxa", how="outer") \
               .merge(df6, on="Taxa", how="outer")


print(merged_df)

merged_df.to_csv("Nitrate_locus_per_taxa.txt", sep="\t", index=False)

df1=pd.read_csv("Nitrate_locus_per_taxa.txt", sep="\t")
df2=pd.read_csv("Taxa_Genero_PW.txt", sep="\t")


# Asumimos que la columna 1 de df1 se llama 'Taxa' y que df2 también tiene una columna 'Taxa'
# Si no tienen nombre, puedes usar df1.columns[0] para referenciarla

# Filtrar los taxones que están en df1 pero no en df2
taxa_df1 = df1[df1.columns[0]]  # primera columna de df1
taxa_df2 = df2[df2.columns[0]]  # primera columna de df2

# Crear máscara de coincidencia
mask = taxa_df1.isin(taxa_df2)

# Filtrar df1 conservando solo los taxones que están en df2
df1_filtered = df1[mask]

# Obtener los taxones eliminados
removed_taxa = df1[~mask][df1.columns[0]].tolist()

# Imprimir los taxones eliminados
print("Taxones eliminados:")
for taxon in removed_taxa:
    print(taxon)

# Guardar el nuevo DataFrame filtrado
df1_filtered.to_csv("Nitrate_locus_filtered.txt", sep="\t", index=False

##############################SULFATE##################################

df1=pd.read_csv("sat_per_specie.txt", sep="\t")
df2=pd.read_csv("aprA_per_specie.txt", sep="\t")
df3=pd.read_csv("aprB_per_specie.txt", sep="\t")
df4=pd.read_csv("dsrA_per_specie.txt", sep="\t")
df5=pd.read_csv("dsrB_per_specie.txt", sep="\t")

merged_df = df1.merge(df2, on="Taxa", how="outer") \
               .merge(df3, on="Taxa", how="outer") \
               .merge(df4, on="Taxa", how="outer") \
               .merge(df5, on="Taxa", how="outer") 


print(merged_df)

merged_df.to_csv("Sulfate_locus_per_taxa.txt", sep="\t", index=False)


df1=pd.read_csv("Sulfate_locus_per_taxa.txt", sep="\t")
df2=pd.read_csv("Taxa_Genero_PW.txt", sep="\t")

# Asumimos que la columna 1 de df1 se llama 'Taxa' y que df2 también tiene una columna 'Taxa'
# Si no tienen nombre, puedes usar df1.columns[0] para referenciarla

# Filtrar los taxones que están en df1 pero no en df2
taxa_df1 = df1[df1.columns[0]]  # primera columna de df1
taxa_df2 = df2[df2.columns[0]]  # primera columna de df2
# Crear máscara de coincidencia
mask = taxa_df1.isin(taxa_df2)
# Filtrar df1 conservando solo los taxones que están en df2
df1_filtered = df1[mask]
# Obtener los taxones eliminados
removed_taxa = df1[~mask][df1.columns[0]].tolist()

# Imprimir los taxones eliminados
print("Taxones eliminados:")
for taxon in removed_taxa:
    print(taxon)

# Guardar el nuevo DataFrame filtrado
df1_filtered.to_csv("Sulfate_locus_filtered.txt", sep="\t", index=False)



