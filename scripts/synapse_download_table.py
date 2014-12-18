import synapseclient
import pandas as pd


def main():
    import argparse
    import sys
    
    parser = argparse.ArgumentParser()

    parser.add_argument('--id', type=str)

    parser.add_argument('outfile', nargs='?',
                        type=argparse.FileType('w'),
                        default=sys.stdout)
    
    args = parser.parse_args()

    syn = synapseclient.login(silent=True)

    qry = syn.queryTable("select * from %s" % args.id)
    res = qry.asDataFrame()
    res.to_csv(args.outfile, sep="\t", index=False)

if __name__ == "__main__":
    main()

