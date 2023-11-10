import argparse
import logging
from tqdm import tqdm



##################################################
# ARGS
##################################################
def get_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', help="Input", required=True)
    parser.add_argument('--ofile', '-o', help="Output file.", required=True)
    return parser


##################################################
# MAIN
##################################################
ARGS = None
if __name__ == "__main__":
    ARGS = get_arg_parser().parse_args()

    pairs = set()
    with open(ARGS.input, 'r') as ifile, open(ARGS.ofile, 'w') as ofile:
        for line in tqdm(ifile):
            sID, qID, dist, pv, sh = line.rstrip('\n').split('\t')
            pair = sorted([sID, qID])
            pair = (pair[0], pair[1])
            # same ID -> skip
            if sID == qID:
                continue
            # group ID already set -> skip
            if pair in pairs:
                continue
            else:
                pairs.add(pair)
                ofile.write('{}\t{}\n'.format(sID, qID))