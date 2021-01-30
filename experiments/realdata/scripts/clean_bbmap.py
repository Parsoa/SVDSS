import sys

import pysam

# sometimes bbmap reports the best alignment as secondary alignment (a secondary alignment has less error than primary one)

def get_errors(al, read_len):
    good_bases = al.get_cigar_stats()[0][7]
    bad_bases = read_len - good_bases
    deletions = al.get_cigar_stats()[0][2]
    return bad_bases + deletions

def main():
    sampath = sys.argv[1]
    sam = pysam.AlignmentFile(sampath, 'r')

    alidx2errors = {}
    read_lens = {} # bbmap doesn't report read sequence on secondary alignment
    for al in sam.fetch():
        read_len = al.query_length
        if al.is_secondary:
            read_len = read_lens[al.query_name]
        else:
            read_lens[al.query_name] = read_len
        errors = get_errors(al, read_len)
        alidx2errors[al.query_name] = errors if al.query_name not in alidx2errors else min(errors,alidx2errors[al.query_name])
    sam.close()

    sam = pysam.AlignmentFile(sampath, 'r')
    print(sam.header, end='')
    for al in sam.fetch():
        read_len = al.query_length
        if al.is_secondary:
            read_len = read_lens[al.query_name]
        else:
            read_lens[al.query_name] = read_len
        errors = get_errors(al, read_len)
        if not al.is_secondary or errors == alidx2errors[al.query_name]:
            print(al.tostring(sam))

if __name__ == "__main__":
    main()
