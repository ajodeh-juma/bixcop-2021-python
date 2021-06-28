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

**Sample Dataset**
```AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC```

**Sample Output**.  
```20 12 17 21```

</details>

<details>
  <summary>Transcribing DNA to RNA</summary>

**Problem**.  
An RNA string is a string formed from the alphabet containing 
```A```, ```C```, ```G```, and ```U```.
Given a DNA string _t_ corresponding to a coding strand, 
its transcribed RNA string _u_ is formed by replacing all 
occurrences of ```T``` in _t_ with ```U``` in _u_.

Given: A DNA string _t_ having length at most 1000 nt.
Return: The transcribed RNA string of _t_

**Sample Dataset**.  
```GATGGAACTTGACTACGTAAATT```

**Sample Output**.  
```GAUGGAACUUGACUACGUAAAUU```
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
If unable to access the file, get it from [here](https://github.com/ajodeh-juma/bixcop-2021-python/blob/main/modules/python/data/test/chr01.fsa)
- Write function(s) to compute the base frequencies.

**Sample Output**.  
Report your output as:<br>
```A: 0.1```<br>```C: 0.2```<br>```G: 0.2```<br>```T: 0.5```
</details>

<details>
  <summary>Parsing GenBank format</summary>
  
**Problem**.  
GenBank format (GenBank Flat File Format) consists of an 
annotation section and a sequence section. The start of the annotation 
section is marked by a line beginning with the word ```LOCUS```. 
The start of sequence section is marked by a line beginning with the 
word ```ORIGIN``` and the end of the section is marked by a 
line with only ```//```.  
Explore a sample genbank file [here](https://www.ncbi.nlm.nih.gov/nuccore/X81322)

**Instructions**.  
Fetch the ```argonaut.gbk``` file []()

**Tasks**.  
1. Write a Python script that parses the GenBank file and performs:  
    * computes sequence records lengths
    * computes GC content
    * reports statistics in an ordered table
2. Use some functionality from the ```BioPython``` package to retrieve the records from GenBank in GenBank format.  
    * Retrieve records for the accessions given in the file [ebov_accessions.txt]() using BioPython Entrez module
    

**Output(s)**.  
- Print a tab-delimited table of accession number, organism name, %GC content, sequence length
- Print the label and sequence of the shortest sequence in FASTA format

**Sample Output**.  
```NM_179453    Arabidopsis thaliana    45.54   3507```.  
```NM_001130718 Strongylocentrotus purpuratus   52.96   2868```.  
```>NM_166020  Drosophila melanogaster ACAGTGCGGAGTGTTTGTTACATGTTAGAGCGTATATATATTTTGAAAAGAGCAGCGACGCCGCCTCAAACCACCGACTAAAATGTCCACGGAGCGTGAGCT```

**Hints**
- Use the Entrez.efetch() function to retrieve the sequences in GenBank format (database “nucleotides”)
- As alternative for your own parser you can also experiment with the Bio.SeqIO.parse() function
- Search field descriptions for sequence database: http://www.ncbi.nlm.nih.gov/books/NBK49540/

**Example**.  
```>>>from Bio import Entrez
>>>Entrez.email = "your_name@your_mail_server.com" 
>>>handle = Entrez.efetch(db="nucleotide", id=["FJ817486, JX069768, JX469983"], rettype="fasta") 
>>>records = handle.read()
>>>print records
```
</details>


 

