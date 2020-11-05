import sys

import pysam

def main():
    sam1_path = sys.argv[1]
    sam2_path = sys.argv[2]

    als = set()
    sam1 = pysam.AlignmentFile(sam1_path, 'r')
    for al in sam1.fetch():
        if al.is_secondary or al.is_unmapped or al.is_supplementary:
            continue
        als.add(al.query_name)

    sam2 = pysam.AlignmentFile(sam2_path, 'r')
    for al in sam2.fetch():
        if al.is_secondary or al.is_unmapped or al.is_supplementary:
            continue
        if al.query_name in als:
            print(al.tostring(sam2))
    

if __name__ == "__main__":
    main()
