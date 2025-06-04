import pandas as pd

def completar_regiones_con_links(df, protein_lengths):
    """
    Completes domain regions by adding linker regions between domains,
    N-terminus, and C-terminus for each protein.
    """
    proteins_dict = {}
    grouped = df.groupby('acc')
    for prot, group in grouped:
        total_length = protein_lengths[prot]
        regions = []
        sorted_features = group.sort_values('dom_start').reset_index(drop=True)
        last_end = 0  # 0-based
        last_feature = None
        for _, row in sorted_features.iterrows():
            start = row['dom_start'] - 1
            end = row['dom_end'] - 1
            current_feature = row['tednum']
            if start > last_end:
                prev_feature = last_feature if last_feature is not None else 'N-terminus'
                next_feature = current_feature
                link_name = (f"{prev_feature}_link_{next_feature}"
                             if last_feature else f"N-terminus_link_{next_feature}")
                regions.append({
                    'acc': prot,
                    'feature': link_name,
                    'start': last_end + 1,
                    'end': start
                })
            regions.append({
                'acc': prot,
                'feature': current_feature,
                'start': start + 1,
                'end': end + 1
            })
            last_end = end + 1
            last_feature = current_feature
        if last_end < total_length:
            link_name = (f"{last_feature}_link_C-terminus"
                         if last_feature else "N-terminus_link_C-terminus")
            regions.append({
                'acc': prot,
                'feature': link_name,
                'start': last_end + 1,
                'end': total_length
            })
        proteins_dict[prot] = regions
    return proteins_dict

def calc_meanf(feature):
    """Return mean of '1' in a binary feature string. If empty, return -1."""
    if feature:
        return feature.count('1') / len(feature)
    else:
        return -1

def get_region_type(string):
    """Classifies region type based on feature name."""
    if 'N-terminus' in string:
        return 'Nter'
    elif 'C-terminus' in string:
        return 'Cter'
    elif 'link' in string:
        d1, _, d2 = string.split("_")
        if d1 == d2:
            return "internal_link"
        else:
            return "link"
    else:
        return 'domain'

def main():
    # File paths
    uniprot_path = '../data/swiss_prot.tsv.gz'
    features_path = '../data/ted_clean/mobidb_features_sequence.tsv'
    ted_domains_path = '../data/ted_clean/ted_domains_sp.tsv'
    output_path = '../data/ted_clean/regions_with_dc_features.tsv'

    print("Loading UniProt table...")
    uniprot = pd.read_table(uniprot_path)
    uniprot.columns = ['acc', 'reviewed', 'name', 'protein_names', 'organism', 'length']
    protein_lengths = uniprot.set_index('acc')['length'].to_dict()

    print("Loading features table...")
    features = pd.read_table(features_path, dtype=str)
    features['length'] = features.acc.map(protein_lengths)

    print("Filling missing values in disorder and missing residues features...")
    fill_defaults = [
        ('prediction-disorder-th_50_sequence', '0'),
        ('prediction-disorder-alphafold_sequence', '0'),
        ('derived-missing_residues-th_90_sequence', '-')
    ]
    for col, fillchar in fill_defaults:
        mask = features[col].isna()
        features.loc[mask, col] = features.loc[mask, 'length'].map(lambda x: fillchar * int(x))

    print("Loading TED domain regions table...")
    ted_domain_regions = pd.read_table(ted_domains_path)
    # Only keep proteins present in features table
    ted_domain_regions = ted_domain_regions[ted_domain_regions.acc.isin(features.acc)]

    print("Building complete protein regions (with linkers)...")
    new_regions = completar_regiones_con_links(ted_domain_regions, protein_lengths)
    protein_regions = pd.concat([pd.DataFrame(i) for i in new_regions.values()]).reset_index(drop=True)

    features_dict = features.set_index('acc').to_dict('index')

    print("Extracting per-region features from sequence annotations...")
    region_feature_cols = [
        ('prediction-disorder-mobidb_lite_sequence', 'lite'),
        ('prediction-disorder-th_50_sequence', 'th50'),
        ('prediction-disorder-alphafold_sequence', 'afd'),
        ('prediction-low_complexity-mobidb_lite_sub_sequence', 'lcd'),
        ('derived-missing_residues-th_90_sequence', 'pdbd')
    ]
    for col, outcol in region_feature_cols:
        protein_regions[outcol] = protein_regions.apply(
            lambda row: features_dict[row.acc][col][row.start - 1: row.end], axis=1
        )

    print("Calculating disorder/LCD/missing residues fractions for each region...")
    for outcol in [f[1] for f in region_feature_cols]:
        dc_col = f"{outcol}_dc"
        protein_regions[dc_col] = protein_regions[outcol].map(calc_meanf)

    protein_regions['length'] = protein_regions['end'] - protein_regions['start'] + 1

    print("Cleaning and filtering invalid regions...")
    # Remove intermediate sequence columns (the original sequences)
    regions = protein_regions.drop(columns=[f[1] for f in region_feature_cols])
    regions = regions.query('th50_dc >= 0').copy()

    print("Assigning region type...")
    regions['tipo'] = regions.feature.map(get_region_type)

    print("Region type summary:")
    print(regions.tipo.value_counts())

    print(f"Saving result to {output_path} ...")
    regions.to_csv(output_path, index=False, sep='\t')
    print("Done.")

if __name__ == '__main__':
    main()
