import pandas as pd

def get_qscat(row):
    """
    Assign quality category based on CATH label granularity and qscore.
    - 'A' if CATH has 4 levels (3 dots) and qscore >= 0.75
    - 'B' if CATH has >=3 levels (>=2 dots) and qscore >= 0.7
    - 'Z' otherwise
    """
    if pd.isnull(row.cath_label) or pd.isnull(row.qscore):
        return "Z"
    if (row.cath_label.count('.') == 3) and (row.qscore >= 0.75):
        return "A"
    elif (row.cath_label.count('.') >= 2) and (row.qscore >= 0.7):
        return "B"
    else:
        return "Z"

def main():
    # File paths
    ted_regions_path = '../data/ted_clean/proteins_ted_regions_df.tsv'
    qscores_path = '../data/ted_ids_qscores.tsv.gz'
    ted_domain_summary_path = '../data/ted_clean/ted_domain_summary_sp.tsv'
    output_path = '../data/ted_clean/ted_domain_summary_sp.tsv'

    print("Loading TED regions...")
    ted_regions = pd.read_table(ted_regions_path)

    print("Loading qscore table...")
    qscores = pd.read_csv(qscores_path, sep='\t', names=['tid', 'qscore'])

    print("Loading TED domain summary...")
    ted_domain_summary = pd.read_table(ted_domain_summary_path)

    print("Merging TED regions with qscore...")
    ted_domain_regions = ted_regions.merge(qscores, on='tid')

    print("Merging with CATH label information...")
    ted_domain_regions = ted_domain_regions.merge(
        ted_domain_summary[['ted_id', 'cath_label']].rename(columns={'ted_id': 'tid'}),
        on='tid',
        how='left'
    )

    print("Assigning quality categories (qcat)...")
    ted_domain_regions['qcat'] = ted_domain_regions.apply(get_qscat, axis=1)

    print(f"Saving results to {output_path} ...")
    ted_domain_regions.to_csv(output_path, sep='\t', index=False)
    print("Done.")

if __name__ == '__main__':
    main()
