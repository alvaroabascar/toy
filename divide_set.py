#!/usr/bin/env python3

import random, sys, argparse

def sample(lista, n):
    """
    Given a list and an integer n, divide this list in two random samples,
    the first of them of size n. Return them as a tuple of two lists.

    """
    # select n indexes, randomly
    first_sample_indexes = random.sample(range(len(lista)), n)
    sample1 = []
    sample2 = []
    for i in range(len(lista)):
        if i in first_sample_indexes:
            sample1.append(lista[i])
        else:
            sample2.append(lista[i])
    return (sample1, sample2)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Take a file and divide it randomly in two new files, the first containing a random sample of the specified size, and the second one containing the rest of the original file.')
    parser.add_argument('-i', '--input',
                        dest = 'input',
                        action = 'store',
                        default = None,
                        help = 'file to split')
    parser.add_argument('-a', '--first-sample',
                        dest = 'first',
                        action = 'store',
                        default = None,
                        help = 'file to store first sample')
    parser.add_argument('-b', '--second-sample',
                        dest = 'second',
                        action = 'store',
                        default = None,
                        help = 'file to store the rest of the original file')

    parser.add_argument('-n', '--first-size',
                        dest = 'size',
                        action = 'store',
                        default = None,
                        help = 'size of the first sample')

    options = parser.parse_args()
    for arg in ['input', 'first', 'second', 'size']:
        if not getattr(options, arg):
            parser.print_usage()
            exit(0)
    with open(options.input, 'r') as stream:
        all_seqs = stream.readlines()
    
    sample1, sample2 = sample(all_seqs, int(options.size))
    
    with open(options.first, 'w') as stream:
        stream.write(''.join(sample1))
    
    with open(options.second, 'w') as stream:
        stream.write(''.join(sample2))
