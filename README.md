# Fast library for finding exact repeat statistics

## **Installation:**

```
git clone https://github.com/at-cg/DiploidGenomeAssembly.git
cd DiploidGenomeAssembly
git clone https://github.com/y-256/libdivsufsort.git 
cd libdivsufsort
mkdir build
cd build && make
echo ‘export LD_LIBRARY_PATH=../libdivsufsort/build/lib:$LD_LIBRARY_PATH’ >> ~/.bashrc
source ~/.bashrc
cd ../../src &&  make
```

## **Test run:**

### Data preparation:
```
python3 prepare.py -mpath <PATH_TO_MATERNAL_HAPLOTYPE_IN_FASTA_FORMAT> -ppath <PATH_WHERE_PATERNAL_HAPLOTYPE_WILL_BE_SAVED_IN_FASTA_FORMAT> -err <SUBSTITUTION_ERROR_RATE>
```

### Computing repeat statistics
```
./summary <PATH_TO_MATERNAL_HAPLOTYPE_IN_FASTA_FORMAT> <PATH_TO_PATERNAL_HAPLOTYPE_IN_FASTA_FORMAT> <PATH_TO_THE_DIRECTORY_THAT_WILL_SAVE_THE_RESULTS>
```


