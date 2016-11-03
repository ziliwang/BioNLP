#!/usr/bin/env python3
import tarfile
import argparse
import os

FUNCTION = '''
batch extract tar.gz
'''
INFO = '''Copyright wzlnot@gmail.com All Rights Reserved. \
Licensed under the MIT License'''


def parse_args():
    parser = argparse.ArgumentParser(
                 description=FUNCTION, epilog=INFO,
                 formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-p', '--path', type=str, required=True,
                              help='path in which extract all *.tar.gz file')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    arg = parse_args()
    paths = []
    for root, dirs, files in os.walk(arg.path):
        for filename in files:
            if os.path.splitext(filename)[-1] in [".gz"]:
                paths.append((os.path.join(root, filename),
                             os.path.join(root, filename.split('.')[0]))
                             )
    for path, target_path in paths:
        tar = tarfile.open(path, "r:gz")
        if not os.path.exists(target_path):
            os.mkdir(target_path)
        tar.extractall(target_path + '/')
        tar.close()
