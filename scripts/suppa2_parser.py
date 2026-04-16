import re
import pandas as pd


def parse_SUPPA2(df, event_type):

    def generate_ID(row):
        try:
            Coord = row["event_id"]

            # split strand
            main_part, strand = Coord.rsplit(":", 1)

            if event_type == "SE":
                Gene, rest = main_part.split(";", 1)
                etype, chrom, up_ee, exon_start, exon_end, down_es = re.split(r"[:-]", rest)

                return f"{Gene};SE,{chrom},{strand};{up_ee}:{exon_start}-{exon_end}:{down_es}"

            elif event_type == "AF":
                Gene, rest = main_part.split(";", 1)
                etype, chrom, c1, c2, c3, c4, c5, c6 = re.split(r"[;:-]", rest)

                if strand == "+":
                    return f"{Gene};AF,{chrom},{strand};{c1}-{c2}:{c3}:{c4}-{c5}:{c6}"
                else:
                    return f"{Gene};AF,{chrom},{strand};{c1}:{c2}-{c3}:{c4}:{c5}-{c6}"

            elif event_type == "AL":
                Gene, rest = main_part.split(";", 1)
                etype, chrom, c1, c2, c3, c4, c5, c6 = re.split(r"[;:-]", rest)

                if strand == "-":
                    return f"{Gene};AL,{chrom},{strand};{c1}-{c2}:{c3}:{c4}-{c5}:{c6}"
                else:
                    return f"{Gene};AL,{chrom},{strand};{c1}:{c2}-{c3}:{c4}:{c5}-{c6}"

            elif event_type == "A5SS":
                Gene, rest = main_part.split(";", 1)
                etype, chrom, c1, c2, c3, c4 = re.split(r"[;:-]", rest)

                if strand == "+":
                    return f"{Gene};A5SS,{chrom},{strand};{c3}-{c1}:{c4}"
                else:
                    return f"{Gene};A5SS,{chrom},{strand};{c1}:{c2}-{c4}"

            elif event_type == "A3SS":
                Gene, rest = main_part.split(";", 1)
                etype, chrom, c1, c2, c3, c4 = re.split(r"[;:-]", rest)

                if strand == "-":
                    return f"{Gene};A3SS,{chrom},{strand};{c3}-{c1}:{c4}"
                else:
                    return f"{Gene};A3SS,{chrom},{strand};{c1}:{c2}-{c4}"

            elif event_type == "RI":
                Gene, rest = main_part.split(";", 1)
                etype, chrom, c1, c2, c3, c4 = re.split(r"[;:-]", rest)

                return f"{Gene};RI,{chrom},{strand};{c1}:{c2}-{c3}:{c4}:{c2}-{c3}"

            elif event_type == "MX":
                Gene, rest = main_part.split(";", 1)
                etype, chrom, c1, c2, c3, c4, c5, c6, c7, c8 = re.split(r"[;:-]", rest)

                return f"{Gene};MX,{chrom},{strand};{c1}:{c2}-{c3}:{c6}-{c7}:{c8}"

            else:
                raise ValueError(f"Unsupported event_type: {event_type}")

        except Exception as e:
            raise ValueError(
                f"Failed to parse SUPPA2 event_id: {row.get('event_id')} | error: {e}"
            )

    # generate uniform_ID
    df = df.copy()
    df["uniform_ID"] = df.apply(generate_ID, axis=1)

    if df["uniform_ID"].isna().any():
        raise ValueError("Some rows failed to generate uniform_ID")

    return df
