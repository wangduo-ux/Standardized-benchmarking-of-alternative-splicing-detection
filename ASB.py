# ASB.py
import os
import re
import argparse
import pandas as pd


from scripts import (
    load_file,
    save_dataframe,
    annotate_dse_class,
    extract_dpsi,
    filter_novel_events,
    plot_venn_per_event_vennlib,
    plot_upset_per_event_combined,
    plot_multiple_event_summaries,
    plot_all_events_grid,
    get_support_software_list,
    ase_extract,
    PARSERS
)


      
def unify_results(input_dir, software_list, event_list, sample_name, output_dir, gtf_path):
    """
    Iterate over each software and event type, parse matching files, and concatenate all results.
    """
    for file_type in ["exclusion_novel", "inclusion_novel"]:
        DSE_number_dict = {}
        DSE_list_dict = {}
        event_dpsi_dict = {}
        test_ase_dict = {}
        control_ase_dict = {}

        for ev in event_list:
            if ev not in ["SE", "A3SS", "A5SS", "MX", "RI", "AF", "AL"]:
                print(f"Warning: unsupported event type '{ev}', skipping.")
                continue

            DSE_list_dict[ev] = {}
            event_dpsi_dict[ev] = {}
            test_ase_dict[ev] = {}
            control_ase_dict[ev] = {}
            summary = []

            support_software = get_support_software_list(ev)
            query_software = [sw2 for sw2 in software_list if any(sw1.lower() in sw2.lower() for sw1 in support_software)]

            for sw in query_software:
                parser_sw = None
                for key, parser in PARSERS.items():
                    if re.search(re.escape(key), sw, re.IGNORECASE):
                        parser_sw = parser
                if parser_sw is None:
                    print(f"Warning: No parser selected for {sw}, skipping.")
                    continue

                df = load_file(input_dir, sample_name, ev, sw)

                if 'psi-sigma' in sw.lower():
                    df = parser_sw(df, ev, gtf_path)
                else:
                    df = parser_sw(df, ev)

                if file_type == "exclusion_novel":
                    df = filter_novel_events(df, sw, ev, input_dir, sample_name)

                if file_type == "inclusion_novel":
                    save_dataframe(df, output_dir, sample_name, sw, ev)

                test_ase_dict[ev][sw], control_ase_dict[ev][sw] = ase_extract(df, sw, ev, sample_name, input_dir)
                df = annotate_dse_class(df, ev, sw)

                df_dpsi = extract_dpsi(df, sw, ev)
                event_dpsi_dict[ev][sw] = df_dpsi[['uniform_ID', 'value']].rename(columns={'value': sw})

                up = len(set(df[df['class'] == 'up-regulate']["uniform_ID"]))
                down = len(set(df[df['class'] == 'down-regulate']["uniform_ID"]))
                summary.append({'Software': sw, 'up-regulate': up, 'down-regulate': down})

                filtered = df[df['class'].isin(['up-regulate', 'down-regulate'])]
                id_set = set(filtered['uniform_ID'])
                if id_set:
                    DSE_list_dict[ev][sw] = id_set

            DSE_number_dict[ev] = pd.DataFrame(summary)

        for data_type, plot_dict in zip(["Test_ASE", "Control_ASE", "DSE"],
                                        [test_ase_dict, control_ase_dict, DSE_list_dict]):
            plot_venn_per_event_vennlib(plot_dict, output_dir, sample_name, filename=f"{data_type}_venn_{file_type}.png")
            plot_upset_per_event_combined(plot_dict, output_dir, sample_name, filename=f"{data_type}_upsets_{file_type}.png")

        plot_multiple_event_summaries(DSE_number_dict, output_dir, sample_name, filename=f"DSE_number_{file_type}.png")
        plot_all_events_grid(event_dpsi_dict, query_software, output_dir, sample_name,
                             filename=f'dpsi_corr_{file_type}', ncols=2, figsize_scale=8)




def main():
    parser = argparse.ArgumentParser(
        description="Unify splicing event results from multiple tools"
    )
    parser.add_argument('--software', '-s',
                        nargs='+',
                        #choices=PARSERS.keys(),
                        required=True,
                        help='space separated list of tools for benchmarking from the following list:SUPPA2 rMATS PSI-Sigma MAJIQ')
    parser.add_argument("-e", "--event", nargs='+', required=False, choices=["SE", "A3SS", "A5SS","MX", "RI", "AF", "AL"],
                        help="list of events to analyze. "
                        "(space separated)\n\n"
                        "Options:\n"
                        "\tSE -- Skipping Exon\n"
                        "\tA3SS -- Alternative Splice Site (3')\n"
                        "\tA5SS -- Alternative Splice Site (5')\n"
                        "\tMX -- Mutually Exclusive Exon\n"
                        "\tRI -- Retained Intron\n"
                        "\tAF -- Alternative First Exon\n"
                        "\tAL -- Alternative Last Exon\n")
    parser.add_argument('--input', '-i',
                        required=True,
                        help='Input file or directory')
    parser.add_argument('--gtf', '-g',
                        required=True,
                        help='a GTF format file')
    parser.add_argument('--output', '-o',
                        required=True,
                        help='Output file')
    parser.add_argument('--sample_name', '-sn',
                        required=True,
                        help='The sample name needs to correspond to the folder name')
    args = parser.parse_args()
    unify_results(args.input, args.software, args.event, args.sample_name, args.output, args.gtf)


if __name__ == '__main__':
    main()
