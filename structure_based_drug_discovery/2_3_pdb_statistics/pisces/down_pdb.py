#!/usr/bin/env python
import sys
import os
from urllib import request
# import gzip
import time

usage = '''
down_pdb.py list_file pdb_dir
'''


def main():
    if len(sys.argv) < 2:
        print(usage)
        sys.exit()

    list_file = sys.argv[1]
    pdb_dir = 'pdb'
    if len(sys.argv) >= 3:
        pdb_dir = sys.argv[2]

    fp = open(list_file)
    lines = fp.readlines()
    fp.close()

    if not os.path.exists(pdb_dir):
        os.makedirs(pdb_dir)

    for line in lines:

        t1 = time.time()
        lis = line.strip().split(' ')
        pdb_id = lis[0]
#        method = lis[1]
#        resolution = lis[2]
#        chain_positions = lis[3]
        pdb_id2 = pdb_id[0:1]
        pdb_dir2 = '%s/%s' % (pdb_dir, pdb_id2)
        if not os.path.exists(pdb_dir2):
            os.makedirs(pdb_dir2)

        line_pdb = 'https://files.rcsb.org/download/%s.pdb.gz' % pdb_id
        pdb_gz_file = '%s/%s.pdb.gz' % (pdb_dir2, pdb_id)
#        pdb_file = '%s/%s.pdb' % (pdb_dir2, pdb_id)

        try:
            request.urlretrieve(line_pdb, pdb_gz_file)
        except Exception as e:
            print(pdb_id, 'error', e)

        t2 = time.time()
        dt = t2 - t1

        print(pdb_id, dt)


if __name__ == "__main__":
    main()
