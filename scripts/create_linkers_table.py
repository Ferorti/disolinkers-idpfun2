import pandas as pd

def load_tables():
    """
    Loads all required dataframes for the analysis.
    """
    print("Loading input tables...")
    arcs = pd.read_table('../data/ted_clean/protein_architectures.tsv')
    features = pd.read_table('../data/ted_clean/mobidb_features_sequence.tsv')
    
    regions = pd.read_table('../data/ted_clean/regions_with_dc_features.tsv')
    regions['tid'] = regions.apply(lambda x: f"AF-{x.acc}-F1-model_v4_{x.feature}", axis=1)
    regions['tidu'] = regions.tid + "_" + regions.start.astype(str)
    
    domains = pd.read_table('../data/ted_clean/ted_domains_sp.tsv')
    domains = domains[domains.acc.isin(features.acc)]

    
    return arcs, features, regions, domains


def proteins_missing_in_regions(domains, regions):
    """
    Prints the number of proteins in domains not present in regions.
    """
    missing = domains[~domains.tid.isin(regions.tid)].drop_duplicates('acc')
    print("Proteins not present in regions dataframe because they don't have mobidb disorder features:")
    print(len(missing))
    return missing

def domains_in_features_not_in_regions(domains, features, regions):
    """
    Returns domains present in features but not in regions.
    """
    dom = domains[domains.acc.isin(features.acc) & ~(domains.tid.isin(regions.tid))]
    print("Domains present in features but missing in regions:", len(dom))
    return dom



def filter_regions(regions, length=10):
    """
    Keeps only valid regions (domain OR length >= length).
    """
    print(f"Filtering regions to keep only 'domain' or length >= {length}...")
    filtered = regions[(regions.tipo == "domain") | (regions.length >= length)].copy()
    return filtered


def merge_architectures(regions, arcs):
    """
    Merges regions with the architectures dataframe.
    """
    print("Merging regions with protein architectures...")
    merged = regions.merge(arcs, on="acc", how="left")
    print("Merged regions with protein architectures.")
    return merged

def extract_linkers(regions):
    """
    Returns DataFrame of linker regions.
    """
    print("Extracting linker regions...")
    return regions[regions.tipo.str.contains('link')].copy()

def main():
    output_path='../data/ted_clean/linkers_features_arcs.tsv'
    
    arcs, features, regions, domains = load_tables()
    
    proteins_missing_in_regions(domains, regions)
    regions = filter_regions(regions)
    regions = merge_architectures(regions, arcs)
    linkers = extract_linkers(regions)
    linkers.to_csv(output_path, sep='\t', index=False)


if __name__ == '__main__':
    main()
