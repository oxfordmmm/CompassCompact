#!/usr/bin/env python
'''
Generate all 0 mask (no mask) file for reference
Usage: python ${COMPASS_ROOT}/nf_ref_nomask.py ${ref} 0
'''

import sys
import pathlib
import argparse

def main(ref_filename,ref_mask_char):
    ref_content = open(ref_filename, 'r').read()
    ref_lines = ref_content.split('\n')

    ref_name = ref_lines[0]
    ref_line_len = len(ref_lines[1])

    ref_base = str(pathlib.Path(ref_filename).stem)
    ref_newmask_filename = ref_base + "_repmask.array"

    with open(ref_newmask_filename, 'w') as f:
        f.write(ref_name + '\n')
        for i in range(2, len(ref_lines)):
            f.write((ref_mask_char * ref_line_len) + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='mask reference with 1 or 0')

    parser.add_argument("-r", dest="ref", required=True,  help="reference Fasta file")
    parser.add_argument("-m", dest="mask", default="0",  help="1 for mask, 0 for no, default 0", choices=['0','1'])

    options = parser.parse_args()

    main(options.ref, options.mask)