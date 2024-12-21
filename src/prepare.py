#!/usr/bin/env python3

import gzip
import io
import random
import click

@click.command()
@click.option('--mat_path', '-mpath', help='path to the maternal FASTA file',
              type=click.Path(exists=True))
@click.option('--pat_path', '-ppath', help='path to the paternal FASTA file',
              type=click.Path(exists=False))
@click.option('--error', '-err', help='substitution error rate', default=1000,
              type=int)
def main(mat_path, pat_path, error):
    s_mat = ''
    with open(mat_path, 'r') as file:
        for line in file:
            if(line[0] == '>'):
                continue
            s_mat += line.strip()

    l = len(s_mat)

    # variant_file_path = '/home/daanish/projects/repeatstatistics/prefix.pair.vcf.gz'

    s_pat = s_mat
    string_list = list(s_pat)
    base = {'A', 'C', 'T', 'G'}
    cnt = 0
    for i in range(l):
        r = random.randint(1, error)
        if(r == 1):
            cnt = cnt + 1
            snp = base.copy()
            snp.remove(string_list[i])
            snplist = list(snp)
            string_list[i] = random.choice(snplist)

    # with gzip.open(variant_file_path, 'rt', encoding='utf-8') as file:
    #     for line in file:
    #         if(line[0] == '#'):
    #             continue
    #         tokens = line.strip().split('\t')
    #         # print(len(tokens[3]), len(tokens[4]), sep = ' ')
    #         if(len(tokens[3]) == 1 and len(tokens[4]) == 1): # SNP
    #             string_list[int(tokens[1]) - 1] = tokens[4][0]

    assert len(s_pat) == len(s_mat)

    # print(cnt) 61197

    s_pat = ''.join(string_list)
    with open(pat_path, 'w') as file:
        file.write('>Modified_string\n')
        file.write(s_pat)


if __name__ == '__main__':
    main()

