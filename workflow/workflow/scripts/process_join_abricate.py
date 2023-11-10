
import pandas as pd

##################################################
# ARGS
##################################################
def get_arg_parser():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', help="Input.", nargs='+', required=True)
    parser.add_argument('--ofile', '-o', help="Output file.", required=True)
    parser.add_argument('--cores', help="Number of cores.", type=int, required=True)
    return parser


##################################################
# MAIN
##################################################
ARGS = None
if __name__ == "__main__":
    ARGS = get_arg_parser().parse_args()

    from os.path import splitext

    dfs = []
    for ifile in ARGS.input:
        d = pd.read_csv(ifile)
        if not d.empty:
            dfs.append(d)
    dfs = pd.concat(dfs, axis=0, sort=False)
    dfs.to_csv(ARGS.ofile, index=False)
    dfs.to_csv(f"{splitext(ARGS.ofile)}.tsv", index=False, sep="\t")
