#!/usr/bin/env python3

import pandas as pd
import os

def process_chunk(chunk):
    """
    Process a chunk of the TED domain boundaries table:
    - Extracts UniProt accession from 'tid'
    - Splits multiple regions in 'positions'
    - Parses domain start/end coordinates
    - Assigns TED domain number
    Returns a DataFrame with one row per domain region.
    """
    try:
        chunk['acc'] = chunk.tid.str.split("-").str[1]
        chunk['positions'] = chunk.positions.str.split('_')
        chunk['tednum'] = chunk.tid.str.split('_').str[-1]
        chunk = chunk.explode('positions')
        chunk['dom_start'] = chunk.positions.str.split('-').str[0].astype(int)
        chunk['dom_end'] = chunk.positions.str.split('-').str[1].astype(int)
        # Filter: only regions where start < end
        chunk = chunk[chunk['dom_start'] < chunk['dom_end']]
    except Exception as e:
        print("Error processing chunk. Exception:", e)
        print(chunk.head())
        raise  # Remove if you want to skip faulty chunks instead of stopping
    return chunk

def main():
    # File paths (edit as needed)
    swissprot_file = '../data/swiss_prot.tsv.gz'
    ted_file = '../data/ted_raw/ted_365m_domain_boundaries_consensus_level.tsv.gz'
    output_dir = '../data/ted_clean'
    output_file = os.path.join(output_dir, 'proteins_ted_regions_df.tsv')

    # Check files exist
    for f in [swissprot_file, ted_file]:
        if not os.path.isfile(f):
            raise FileNotFoundError(f"File not found: {f}")

    # Create output dir if missing
    os.makedirs(output_dir, exist_ok=True)

    # Load Swiss-Prot accessions
    print("Loading Swiss-Prot accessions...")
    uniprot_sp = pd.read_table(
        swissprot_file,
        names=['acc', 'reviewed', 'name', 'protein_names', 'organism', 'length']
    )
    sp_accs = set(uniprot_sp.acc.drop_duplicates())

    # Process TED boundaries in chunks
    print("Processing TED domain boundaries...")
    chunks = pd.read_table(
        ted_file,
        chunksize=1_000_000,
        header=None,
        names=['tid', 'positions', 'level']
    )

    proteins_ted_regions_list = []
    n_chunks = 0
    n_rows = 0

    for chunk in chunks:
        print(f"Processing chunk {n_chunks}, total rows so far: {n_rows}", end='\r')
        chunk = process_chunk(chunk)
        chunk = chunk[chunk.acc.isin(sp_accs)]
        n_chunks += 1
        n_rows += len(chunk)
        proteins_ted_regions_list.append(chunk)

    print(f"\nMerging results ({n_rows} regions in total)...")
    proteins_ted_regions_df = pd.concat(proteins_ted_regions_list, ignore_index=True)

    # Save result
    print(f"Saving output to {output_file} ...")
    proteins_ted_regions_df.to_csv(output_file, sep='\t', index=False)
    print("Done.")

if __name__ == '__main__':
    main()
