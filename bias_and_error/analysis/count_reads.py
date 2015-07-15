# this script counts the number of reads for each transcript based on the
# eXpress format

from collections import defaultdict
import io
import sys

def get_read_info(line):
    line = line.split(":")
    return {
        'rid': line[0],
        'target_id': line[1],
            }

def main():
    # take one read
    fname = sys.argv[1]
    line_num = -1

    hist = defaultdict(int)

    with open(fname, "r") as fhandle:
        for line in fhandle:
            # skip all lines other than the ones containing sim info
            line_num += 1
            if line_num % 2 != 0:
                continue
            read_info = get_read_info(line)
            hist[read_info['target_id']] += 1

    sorted_keys = sorted(hist.keys())
    vals = []
    for k in sorted_keys:
        print k, "\t", hist[k]


if __name__ == '__main__':
    main()
