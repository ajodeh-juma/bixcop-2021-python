# Bioinformatics community of Practice - Python programming

[![Twitter Follow](https://img.shields.io/twitter/follow/john_juma.svg?style=social)](https://twitter.com/john_juma)


[![forthebadge made-with-python](http://ForTheBadge.com/images/badges/made-with-python.svg)](https://www.python.org/)
[![Made withJupyter](https://img.shields.io/badge/Made%20with-Jupyter-orange?style=for-the-badge&logo=Jupyter)](https://jupyter.org/try)

```python``` programming exercises for the Bioinformatics Community of Practice 2021,
hosted at the [International Livestock Research Institute](https://www.ilri.org/), Nairobi, Kenya.    
Largely inspired by [Rosalind](http://rosalind.info/about/) problems.  

**Instructions**.  
1. Download and Install ```anaconda``` installer from the [website](https://www.anaconda.com/products/individual).  
2. Create a ```conda``` environment, and provide it with a simple name e.g ```python-programming-env```.  
3. Install the packages:
    * [biopython](https://biopython.org/)
    * [fastx_toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)

**Requirements**.  
    - ```anaconda/miniconda```
    - ```jupyter```
    - ```biopython```
    - ```fastx_toolkit```
    
**References**.  
1.  Rosalind: http://rosalind.info/about/
2. fastx_toolkit: The FASTX-Toolkit is a collection of command line tools for Short-Reads FASTA/FASTQ files preprocessing. 
3. Cock _et al_. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics.  _Bioinformatics_, Volume **25**, Issue 11, 1 June 2009, Pages 1422–1423, [DOI](https://doi.org/10.1093/bioinformatics/btp163)

![](https://img.shields.io/badge/licence-MIT-lightgrey.svg)

**Tasks**
<details>
  <summary>Variables and arithmetic</summary>
  
**Problem**.  

**Given**: Two positive integers ``a`` and ```b```, each less than 1000.

**Return**: The integer corresponding to the square of the hypotenuse of the right triangle whose legs have lengths ```a```
and ```b```.  

**Sample Dataset**.  
```3 5```

**Sample Output**.  
```34```

##### Instructions
- Obtain the [dataset](https://github.com/ajodeh-juma/bixcop-2021-python/raw/main/data/test/ini2_dataset.txt) 
using ```wget``` command
- Write a function to print the output

</details>

<details>
  <summary>Counting DNA Nucleotides</summary>
  
**Problem**.  
A string is simply an ordered collection of symbols selected from 
some alphabet and formed into a word; the length of a string is the 
number of symbols that it contains.

An example of a length 21 DNA string (whose alphabet contains the 
symbols ```A```, ```C```, ```G```, and ```T```) is 
```ATGCTTCAGAAAGGTCTTACG```.

Given: A DNA string _s_ of length at most 1000 nt.  

Return: Four integers (separated by spaces) counting the 
respective number of times that the symbols ```A```, ```C```, ```G```, 
and ```T``` occur in _s_.

**Sample Dataset**.  
```AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC```

**Sample Output**.  
```20 12 17 21```

##### Instructions
- Obtain the [dataset](https://github.com/ajodeh-juma/bixcop-2021-python/raw/main/data/test/dna_dataset.txt)
using ```wget``` command
- Write a function to print the output
</details>

<details>
  <summary>Transcribing DNA to RNA</summary>

**Problem**.  
An RNA string is a string formed from the alphabet containing 
```A```, ```C```, ```G```, and ```U```.
Given a DNA string _t_ corresponding to a coding strand, 
its transcribed RNA string _u_ is formed by replacing all 
occurrences of ```T``` in _t_ with ```U``` in _u_.

**Given**: A DNA string _t_ having length at most 1000 nt.  

**Return**: The transcribed RNA string of _t_

**Sample Dataset**.  
```GATGGAACTTGACTACGTAAATT```

**Sample Output**.  
```GAUGGAACUUGACUACGUAAAUU```

##### Instructions
- Obtain the [dataset](https://github.com/ajodeh-juma/bixcop-2021-python/raw/main/data/test/rna_dataset.txt),
using ```wget``` command
- Write a function to print the output
</details>

<details>
  <summary>Complementing a Strand of DNA</summary>

**Problem**.  
In DNA strings, symbols ```A``` and ```T``` are complements of each other, as are ```C``` and ```G```.  

The reverse complement of a DNA string _s_ is the string _s_<sup>_c_</sup> formed by reversing the symbols of _s_,
then taking the complement of each symbol (e.g., the reverse complement of ```GTCA``` is ```TGAC```).

**Given**: A DNA string _s_ of length at most 1000 bp.  

**Return**: The reverse complement _s_<sup>_c_</sup> of _s_.

**Sample Dataset**.  
```AAAACCCGGT```

**Sample Output**.  
```ACCGGGTTTT```

##### Instructions
- Obtain the [dataset](https://github.com/ajodeh-juma/bixcop-2021-python/raw/main/data/test/revc_dataset.txt),
using ```wget``` command
- Write a function to print the output
</details>

<details>
  <summary>Counting point mutations</summary>

**Problem**.  
Given two strings _s_ and _t_ of equal length, the Hamming distance between _s_ and _t_, 
denoted _d_<sub>H</sub>(_s_,_t_), is the number of corresponding symbols that differ in _s_ and _t_, see the figure below.  
 
![Figure 2](https://github.com/ajodeh-juma/bixcop-2021-python/raw/main/images/Hamming_distance.png).  

Given: Two DNA strings _s_ and _t_ of equal length (not exceeding 1 kbp).

Return: The Hamming distance _d_<sub>H</sub>(_s_,_t_)

**Sample Dataset**.  
```
GAGCCTACTAACGGGAT
CATCGTAATGACGGCCT
```

**Sample Output**.  
```7```

##### Instructions
- Obtain the [dataset](https://github.com/ajodeh-juma/bixcop-2021-python/raw/main/data/test/hamm_dataset.txt),
using ```wget``` command
- Write a function to print the output
</details>

<details>
  <summary>Computing GC content</summary>

**Problem**.  
The GC-content of a DNA string is given by the percentage of symbols in the string that are ```C``` or ```G```.  
For example, the GC-content of ```AGCTATAG``` is ```37.5%```.  Note that the reverse complement of any DNA string has the same GC-content.

DNA strings must be labeled when they are consolidated into a database.  
A commonly used method of string labeling is called ```FASTA``` format.  
In this format, the string is introduced by a line that begins with ```>```, 
followed by some labeling information. Subsequent lines contain the string itself; the first line to begin with ```>``` 
indicates the label of the next string.

**Given**: At most 10 DNA strings in ```FASTA``` format (of length at most 1 kbp each).

**Return**: The ID of the string having the highest GC-content, followed by the GC-content of that string in 6 decimal places. 


**Sample Dataset**.  
```
>Rosalind_6404
CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
TCCCACTAATAATTCTGAGG
>Rosalind_5959
CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
ATATCCATTTGTCAGCAGACACGC
>Rosalind_0808
CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC
TGGGAACCTGCGGGCAGTAGGTGGAAT
```

**Sample Output**.  
```
Rosalind_0808
60.919540
```

##### Instructions
- Obtain the [dataset](https://github.com/ajodeh-juma/bixcop-2021-python/raw/main/data/test/gc_dataset.txt),
using ```wget``` command
- Write a function to print the output
</details>


<details>
  <summary>Finding Base Frequencies</summary>

**Problem**.  
DNA consists of four molecules called nucleotides, or bases, and can be 
represented as a string of the letters ```A```, ```C```, ```G```, 
and ```T```. But this does not 
mean that all four nucleotides need to be similarly frequent. 
Are some nucleotides more frequent than others, say in yeast, as 
represented by the first chromosome of yeast? Also, DNA is really not a 
single thread, but two threads wound together. This wounding is based on 
an ```A``` from one thread binding to a ```T``` of the other thread, 
and ```C``` binding to 
```G``` (that is, ```A``` will only bind with ```T```, not with ```C``` or ```G```). 
Could this fact 
force groups of the four symbol frequencies to be equal? 
The answer is that the A-T and G-C binding does not in principle force 
certain frequencies to be equal, but in practice they usually become so 
because of evolutionary factors related to this pairing.

##### Task
Compute the frequencies of the bases ```A```, ```C```, ```G```, 
and ```T```. That is, the number of times each base occurs in the 
DNA string, divided by the length of the string. For example, 
if the DNA string is ```ACGGAAA```, the length is 7, ```A``` appears 
4 times with frequency 4/7, ```C``` appears once with frequency 1/7, 
```G``` appears twice with frequency 2/7, and ```T``` does not 
appear so the frequency is 0.

##### Instructions
- Obtain the [yeast chromosome 1 sequence](http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr01.fsa),
If unable to access the file, get it from [here](https://github.com/ajodeh-juma/bixcop-2021-python/raw/main/data/test/chr01.fsa)
- Write function(s) to compute the base frequencies.

**Sample Output**.  
(Note that the numbers in this output are made up. The numbers in your output will be different!).  

```A: 0.112647```<br>```C: 0.223456```<br>```G: 0.212356```<br>```T: 0.501349```.  


</details>

<details>
  <summary>Parsing GenBank format file</summary>
  
**Problem**.  
GenBank format (GenBank Flat File Format) consists of an 
annotation section and a sequence section. The start of the annotation 
section is marked by a line beginning with the word ```LOCUS```. 
The start of sequence section is marked by a line beginning with the 
word ```ORIGIN``` and the end of the section is marked by a 
line with only ```//```.  
Explore a sample genbank file [here](https://www.ncbi.nlm.nih.gov/nuccore/X81322)

**Instructions**.  
Fetch the ```argonaut.gbk``` [here](https://github.com/ajodeh-juma/bixcop-2021-python/raw/main/data/test/argonaut.gb)
using the ```wget``` command.

**Tasks**.  
1. Write your own Python script that parses the GenBank file and performs:  
    * computes sequence records lengths
    * computes GC content
    * reports statistics in an ordered table
2. Use some functionality from the ```BioPython``` package to retrieve the records from GenBank in GenBank format.  
    * Retrieve records for the accessions given in the file [ebov_accessions.txt](https://github.com/ajodeh-juma/bixcop-2021-python/raw/main/data/test/ebov_accessions.txt) 
    using ```BioPython Entrez``` module.  
    * Save your genbank records as ```ebov.gbk```
    
        **Hints**.  
        * Use the ```Entrez.efetch()``` function to retrieve the sequences in GenBank format (database “nucleotides”)
        * As alternative for your own parser you can also experiment with the ```Bio.SeqIO.parse()``` function
        * Search field descriptions for sequence database: http://www.ncbi.nlm.nih.gov/books/NBK49540/.  
        
         **Example**.  
        ```
        >>>from Bio import Entrez
        >>>Entrez.email = "your_name@your_mail_server.com" 
        >>>handle = Entrez.efetch(db="nucleotide", id=["FJ817486, JX069768, JX469983"], rettype="fasta") 
        >>>records = handle.read()
        >>>print(records)
        ```
    

**Output(s)**.  
* Print a tab-delimited table of accession number, organism name, %GC content, sequence length
* Print the label and sequence of the shortest sequence in ```FASTA``` format

**Sample Output**.  
```NM_179453    Arabidopsis thaliana    45.54   3507```.  
```NM_001130718 Strongylocentrotus purpuratus   52.96   2868```.  
```>NM_166020  Drosophila melanogaster ACAGTGCGGAGTGTTTGTTACATGTTAGAGCGTATATATATTTTGAAAAGAGCAGCGACGCCGCCTCAAACCACCGACTAAAATGTCCACGGAGCGTGAGCT```

</details>

<details>
  <summary>Trimming NGS data</summary>
  
  
**Description**.  
Next-generation sequencing machines produce vast amounts of DNA or RNA reads. 
Illumina’s sequencers produce output in the form of FASTQ files. 
Quality control of the produced reads is a necessary step before any downstream analysis, 
such as assembly or mapping. Typically, the average quality at the 3’ end of 
the reads is lower than at the 5’ end of the read, caused by the 
chemistry and process of sequencing. When plotting the ```per-base quality``` 
for all reads in a FASTQ file, a typical pattern looks like this:  

![Figure 2](https://github.com/ajodeh-juma/bixcop-2021-python/raw/main/images/per_base_quality.png).  


In this assignment you will calculate the average quality score at each position of the read, 
use a command line tool to ```trim``` off low-quality bases, and assess whether the average per-base quality has improved.

**Assignment**.  
Write a script that performs the following tasks:
1. Parse a FASTQ file. Translate the quality values (```Illumina 1.5+ encoding```) to a
scale from ```0``` to ```41```. Use the built-in Python function ```ord()``` for the translation.
2. Calculate the length of the shortest sequence in the input FASTQ file, the longest
sequence in the file, and calculate the average sequence length.
3. Calculate the average quality score (on a scale from ```0``` to ```41```) at each position of
the read. In the raw FASTQ file all sequences have the same length. When calculating the average quality value at position 0, you average over the quality values at position 0 of all the reads, etc. Your script should be able to work on input sequences of any length (e.g. the tiny example below).
4. In your Python script, trim off low-quality bases using the program ```fastq_quality_trimmer```. Set the quality threshold to ```30```. Name the output file ‘trimmed.fq’.
5. Calculate the average quality score at each position of the read in the trimmed file. At each position, calculate the improvement with respect to that of the original FASTQ file (step 3).
6. Report the minimum, maximum, average sequence length for both the original FASTQ file, and the trimmed FASTQ file.
7. For each read position, report the average quality score in the original file, the average quality score in the trimmed file, and the improvement in average quality, in tab-delimited columns.
Input
A FASTQ file containing 10000 records. It is a sample of genomic reads from a tomato plant.  
Use the command ```wget``` to download the file [here](https://github.com/ajodeh-juma/bixcop-2021-python/raw/main/data/test/tomatosample.fq).  
For development purposes or if you fail to get step 4 working, 
a trimmed FASTQ file is also provided: http://www.bioinformatics.nl/courses/BIF-30806/docs/trimmedsample.fq

**Output**.  
The output of your script should look like this: 
(Note that the numbers in this output are made up. The numbers in your output will be different!).  
```ORIGINAL: min=100, max=100, avg=100.00```<br>```TRIMMED: min=27, max=100, avg=96.48```.  
```1    33.18   0.00```<br>```2 33.42   0.00```<br>```3 33.11   0.00```<br>```99   28.22   0.00```<br>```100    27.82   0.00```.  

**A tiny example**.  

| Label         | Original      | Quality in original   |  Trimmed at t=30  | Quality in trimmed |
|:-------------:|:-------------:| :-------------------: | :----------------:| :-----------------:|
| Seq1          | AGACA         | 34,34,34,37,37        | AGACA             | 34,34,34,37,37     |
| Qual1         | bbbee         |                       | bbbee             |                    |
| Seq2          | CCCAA         | 40,40,40,39,27        | CCCA              | 40,40,40,39        |
| Qual2         | hhhg[         |                       | hhhg              |                    |
| Seq3          | ATAAT         | 35,35,35,3,2          | ATA               | 35, 35, 35         |
| Qual3         | cccCB         |                       | ccc               |                    |
|               | pos 1 avg     |     36.33             |                   |      36.33         |
|               | pos 4 avg     |     26.33             |                   |      38.00         |


The tiny example used can be obtained for validation purposes: http://www.bioinformatics.nl/courses/BIF-30806/docs/tiny.fq.  
**Environment**.  
- create a ```conda environment``` named ```bioinfm-env```.  
- install the package ```fastx_toolkit``` and use the program ```fastq_quality_trimmer``` to ```trim``` sequences.  
- Try it by typing ```fastq_quality_trimmer –h```.  
- On the command line. You should see information on the usage and options.  

**Additional notes**:  
- Put your full name and student number as a comment in your script and put your name in the file name of your script.
- You may use the slides and the code from your exercises from this week. You cannot use BioPython or comparable packages, you should write your own fastq parser. You may not directly copy code from the internet, but you may use it as inspiration.
- The FASTQ format and quality values are explained on: http://en.wikipedia.org/wiki/FASTQ_format
- Running ```fastq_quality_trimmer``` takes only a few seconds or so. But to avoid running it over and over again, make sure your code checks whether the trimmed file exists.
- Think about your code organization. Use subroutines.
- Document your code.
- Make sure you hand in a working script. If it is unfinished, you can leave the unfinished part in comments.
</details>

