from argparse import ArgumentParser
import sys


# Allen McBride
# June 5, 2024

# Given an iterable container of strings, each representing a residue number,
# construce a "restraintmask" list, in the format required by AmberTools'
# sander tool. This restraintmask designates names those residues NOT in the
# input container. The intent is that resNumberStrings represents all residues
# involved in clashes, and the resulting restraintmask restrains all OTHER
# residues during minimization.

# If invoked from the command line:
#
# Standard input should be a sequence of lines, each containing a string
# representing a residue number.
#
# --maxres <res> (<res> is the highest-numbered residue in the structure to be
# minimized)


def findMask(resNumberStrings, maxRes):
    skiplist = []
    for line in resNumberStrings:
        skiplist.append(int(line))

    waitingForFirst = True
    mask = ':'
    for n in range(1, maxRes + 1):
        if waitingForFirst:
            if n not in skiplist:
                waitingForFirst = False
                first = n
        else:
            if n in skiplist:
                second = n - 1
                if second == first:
                    mask += f'{first},'
                else:
                    mask += f'{first}-{second},'
                waitingForFirst = True
    if not waitingForFirst:
        mask += f'{first}-{maxRes}'
    else:
        mask.rstrip(',')
    return(mask)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("--maxres", required=True, type=int)
    args = parser.parse_args()
    print(findMask(sys.stdin, args.maxres))
