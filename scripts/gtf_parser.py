import pandas as pd

def extract_gene_transcript_strand(gtf_file: str) -> pd.DataFrame:
    """
    Extract gene ID, transcript ID, and strand information from a GTF file.

    Parameters
    ----------
    gtf_file : str
        Path to the input GTF file.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing three columns: gene_id, transcript_id, and strand.
    """
    records = []
    with open(gtf_file, 'r') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            info = fields[8]
            strand = fields[6]
            gene_id = transcript_id = None
            for item in info.strip().split(';'):
                item = item.strip()
                if item.startswith('gene_id'):
                    gene_id = item.split()[1].strip('"')
                elif item.startswith('transcript_id'):
                    transcript_id = item.split()[1].strip('"')
            if gene_id and transcript_id:
                records.append((gene_id, transcript_id, strand))

    df = pd.DataFrame(records, columns=['gene_id', 'transcript_id', 'strand'])
    return df.drop_duplicates()
