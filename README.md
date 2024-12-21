# Fast library for finding exact repeat statistics

## **Dependencies:**
> c++17

## **Installation:**
```
git clone https://github.com/at-cg/DiploidGenomeAssembly.git
cd DiploidGenomeAssembly
git clone https://github.com/y-256/libdivsufsort.git 
cd libdivsufsort
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE="Release" -DCMAKE_INSTALL_PREFIX="/usr/local" ..
make
echo export LD_LIBRARY_PATH="../libdivsufsort/build/lib:$LD_LIBRARY_PATH" >> ~/.bashrc
source ~/.bashrc
cd ../../src &&  make
```

## **Test run:**

### Data preparation:
```
python3 prepare.py -mpath <PATH_TO_MATERNAL_HAPLOTYPE_IN_FASTA_FORMAT> -ppath <PATH_WHERE_PATERNAL_HAPLOTYPE_WILL_BE_SAVED_IN_FASTA_FORMAT> -err <INVERSE_OF_SUBSTITUTION_ERROR_RATE (must be an integer, default = 1000)>
```

### Computing repeat statistics
```
./summary <PATH_TO_MATERNAL_HAPLOTYPE_IN_FASTA_FORMAT> <PATH_TO_PATERNAL_HAPLOTYPE_IN_FASTA_FORMAT> <PATH_TO_THE_DIRECTORY_THAT_WILL_SAVE_THE_RESULTS> <THRESHOLD ON THE MINIMUM REPEAT LENGTH AMONGST ALL POSSIBLE TYPE OF REPEATS (must be an integer)> <WHETHER INTERLEAVED_STATS ARE REQUIRED (0 (no) / 1 (yes))>
```

**NOTE:** Computation of interleaved statistics may take some time depending on the size of the input. The result for the test example below can be reproduced in a few seconds.

## **Example:**
```
python3 prepare.py -mpath ../test/data/chrMT.fna -ppath ../test/data/chrMT_m.fna -err 1000
./summary ../test/data/chrMT.fna ../test/data/chrMT_m.fna ../results 1 1
```

