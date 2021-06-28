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
 

