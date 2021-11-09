from utils  import run_epost_split
import argparse
import sys

parser = argparse.ArgumentParser(description='Query information using eutils')


parser.add_argument('--log_file',type=str, nargs='?', default="pipeline.log")
parser.add_argument('--df_file',type=str)
parser.add_argument('--ofile',type=str)
parser.add_argument('--header',type=str)
parser.add_argument('--cmd',type=str)
parser.add_argument('--df_col',type=str)
parser.add_argument('--split_size', type=int)
parser.add_argument('--split_str',type=str,nargs='?')

parser.add_argument('--db_source', nargs='?',type=str)
parser.add_argument('--db_target', nargs='?',type=str)

args=vars(parser.parse_args())
run_epost_split(**args)