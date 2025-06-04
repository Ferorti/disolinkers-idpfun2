import pandas as pd
import json
import os

def load_mobidb_jsonl(filepath):
    """Load a MobiDB JSONL file into a dict."""
    data = {}
    with open(filepath) as f:
        for line in f:
            json_data = json.loads(line.strip())
            data[json_data.get('acc')] = json_data
    return data

def load_uniprot(filepath):
    """Load Swiss-Prot table into a DataFrame."""
    df = pd.read_table(filepath, header=None)
    df.columns = ['acc', 'reviewed', 'name', 'protein_names', 'organism', 'length']
    return df

def region_to_seq(regions, seq_len, neg_char='0', pos_char='1'):
    """Convert list of (start, end) to binary sequence string."""
    base = [neg_char] * seq_len
    for region in regions:
        for pos in range(region[0], region[1]):
            if pos > seq_len:
                return "error"
            base[pos - 1] = pos_char
    return "".join(base)

def fake_no_regions(length, neg_char='0'):
    """Generate a negative-character-only sequence."""
    return neg_char * length

def extract_seq_features(data, features):
    """Generate sequence features for all proteins in data."""
    seq_features = {}
    for acc, protein_data in data.items():
        new_data = {}
        for feature in features:
            if feature in protein_data:
                regions = protein_data[feature].get('regions')
                if regions:
                    feature_seq = region_to_seq(regions, protein_data['length'])
                else:
                    feature_seq = fake_no_regions(protein_data['length'])
                new_data[feature + "_sequence"] = feature_seq
        seq_features[acc] = new_data
    return seq_features

def fill_missing_features(df, largos):
    """Fill missing values in feature columns."""
    for feature, neg_char in [
        ('prediction-disorder-mobidb_lite_sequence', '0'),
        ('prediction-disorder-th_50_sequence', '0'),
        ('prediction-low_complexity-mobidb_lite_sub_sequence', '0')
    ]:
        df[feature] = df[feature].fillna(df.acc.map(lambda x: neg_char * largos[x]))
    df['derived-missing_residues-th_90_sequence'] = df['derived-missing_residues-th_90_sequence'].fillna(
        df.acc.map(lambda x: '-' * largos[x])
    )
    return df

def main():
    mobidb_jsonl = '../data/reduced_mobidb.mjson'
    swissprot_file = '../data/swiss_prot.tsv.gz'
    output_file = '../data/ted_clean/mobidb_features_sequence.tsv'

    print("Loading MobiDB data...")
    data = load_mobidb_jsonl(mobidb_jsonl)
    print("Loading Swiss-Prot...")
    uniprot = load_uniprot(swissprot_file)

    features = [
        'prediction-disorder-mobidb_lite',
        'prediction-disorder-th_50',
        'derived-missing_residues-th_90',
        'prediction-disorder-alphafold',
        'prediction-plddt-alphafold',
        'prediction-low_complexity-mobidb_lite_sub'
    ]

    print("Extracting sequence features...")
    seq_features = extract_seq_features(data, features)

    print("Generating DataFrame...")
    seq_features_df = pd.DataFrame.from_dict(seq_features, orient='index').reset_index(names='acc')

    largos = {k: v['length'] for k, v in data.items()}
    seq_features_df = fill_missing_features(seq_features_df, largos)

    print(f"Saving table to {output_file}...")
    seq_features_df.to_csv(output_file, sep='\t', index=False)
    print("Done.")

if __name__ == '__main__':
    main()
