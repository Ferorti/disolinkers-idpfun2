import pandas as pd

def preprocess_regions(regions_path, domain_regions_path):
    """
    Loads regions and TED domain regions, constructs 'tid', and merges CATH label and qcat.
    """
    print("Loading region and TED domain tables...")
    regions = pd.read_table(regions_path)
    ted_domain_regions = pd.read_table(domain_regions_path)
    regions['tid'] = regions.apply(lambda x: f'AF-{x.acc}-F1-model_v4_{x.feature}', axis=1)
    merge_cols = ted_domain_regions[['tid', 'cath_label', 'qcat']].drop_duplicates()
    regions = regions.merge(merge_cols, how='left')
    return regions

def categorize_disorder_columns(regions):
    """
    Adds categorical columns for disorder/LCD features based on defined bins and labels.
    """
    print("Categorizing disorder features into bins...")
    bins = [0, 0.7, 0.8, 0.9, 1]
    labels = list("0321")  # Higher = more disordered
    for col, cat_col in [
        ('th50_dc', 'th50_cat'),
        ('lite_dc', 'lite_cat'),
        ('afd_dc', 'afd_cat')
    ]:
        regions[cat_col] = pd.cut(regions[col], bins=bins, labels=labels, include_lowest=True)
    regions['afd_dc_cat'] = regions.afd_dc.fillna("0")
    return regions

def build_cat_arq(df, col="th50_cat"):
    """
    Builds architecture string for a protein, using type/cat columns.
    Converts everything to string, using '-' for missing values.
    """
    arq = []
    for _, row in df.iterrows():
        if row.tipo == 'Nter' or row.tipo == 'Cter':
            a = '_'
        elif "link" in row.tipo:
            a = row[col]
        else:
            a = row.qcat
        if pd.isnull(a):
            a = '-'
        arq.append(str(a))
    return "".join(arq)

def build_protein_architectures(regions, min_lengths, cat_cols):
    """
    For each protein and each minimum region length, builds architecture strings
    using different categorical columns.
    """
    print("Building architecture strings for each protein...")
    results = []
    for i, (acc, group) in enumerate(regions.groupby("acc"), 1):
        if i % 500 == 0:
            print(f"Processed {i} proteins...")
        row_data = {"acc": acc}
        for min_len in min_lengths:
            sub = group.query("length >= @min_len")
            for cat_col in cat_cols:
                suffix = cat_col.replace("_dc_cat", "").replace("_cat", "")
                arq_col_name = f"{suffix}_cat_{min_len:02d}"
                row_data[arq_col_name] = build_cat_arq(sub, col=cat_col)
        results.append(row_data)
    return pd.DataFrame(results)

def main():
    regions_path = '../data/ted_clean/regions_with_dc_features.tsv'
    domain_regions_path = '../data/ted_clean/ted_domains_sp.tsv'
    output_path = '../data/ted_clean/protein_architectures.tsv'
    min_lengths = [10, 20]
    cat_cols = ["th50_cat", "lite_cat"]

    # Load, merge, and categorize
    regions = preprocess_regions(regions_path, domain_regions_path)
    regions = regions.drop(columns='tid')
    regions = categorize_disorder_columns(regions)
    regions = regions.query('th50_dc >= 0').copy()

    # Build architectures
    proteins_arc = build_protein_architectures(regions, min_lengths, cat_cols)

    print(f"Saving result to {output_path} ...")
    proteins_arc.to_csv(output_path, sep='\t', index=False)
    print("Done.")

if __name__ == '__main__':
    main()
