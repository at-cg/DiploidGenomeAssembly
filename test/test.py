import random
import click

@click.command()
@click.option('--num', '-num', default=1, help='number of haplotypes',
              type=int)
@click.option('--len', '-len', help='length of each haplotype',
              type=int)
@click.option('--snp', '-snp', default=10, help='probability with which SNP occur',
              type=int)
@click.option('--output_path1', '-opath1', help='path to the file that will store the fasta file for haplotype1', 
              type=click.Path(exists=False))
@click.option('--output_path2', '-opath2', help='path to the file that will store the fasta file for haplotype2', 
              type=click.Path(exists=False))
def main(num, len, snp, output_path1, output_path2):
    s = ''
    d = {0 : 'A', 1 : 'C', 2 : 'G', 3 : 'T'}

    for i in range(len):
        s += d[random.randint(0, 3)]

    with open(output_path1, 'w') as file:
        print(s, file = file)

    if num == 2:
        s1 = s
        string_list = list(s1)
        base = {'A', 'C', 'T', 'G'}
        for i in range(len):
            r = random.randint(1, snp)
            if(r == 1):
                baselist = base.copy()
                baselist.remove(string_list[i])
                snplist = list(baselist)
                string_list[i] = random.choice(snplist)

        s1 = ''.join(string_list)
        with open(output_path2, 'w') as file:
            file.write(s1)


if __name__ == '__main__':
    main()
