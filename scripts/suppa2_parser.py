import pandas as pd

def parse_SUPPA2(df: pd.DataFrame, event_type: str) -> pd.DataFrame:
    """
    Parse SUPPA2 event results based on the specified event type.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame containing SUPPA2 results, with event IDs in the index.
    event_type : str
        Event type to filter. Supported types include:
        - SE  : Skipped Exon
        - A3SS: Alternative 3' Splice Site
        - A5SS: Alternative 5' Splice Site
        - AF  : Alternative First Exon
        - AL  : Alternative Last Exon
        - RI  : Retained Intron
        - MX  : Mutually Exclusive Exon

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame with a new column `uniform_ID` containing event IDs
        that match the specified event type.

    Raises
    ------
    ValueError
        If an unsupported event_type is provided.
    """
    event_map = {
        "SE": "SE", "A3SS": "A3", "A5SS": "A5",
        "AF": "AF", "AL": "AL", "RI": "RI", "MX": "MX"
    }
    event = event_map.get(event_type)
    if not event:
        raise ValueError(f"Unsupported event_type: {event_type}")

    return (
        df.reset_index()
          .rename(columns={'index': 'uniform_ID'})
          .query("uniform_ID.str.contains(@event)")
    )
