##### Finding Base Frequencies

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

#### Task
Compute the frequencies of the bases ```A```, ```C```, ```G```, and ```T```. 
That is, the number of times each base occurs in the DNA string, 
divided by the length of the string. For example, 
if the DNA string is ```ACGGAAA```, the length is 7, ```A``` appears 
4 times with frequency 4/7, ```C``` appears once with frequency 1/7, 
```G``` appears twice with frequency 2/7, and ```T``` does not appear so the 
frequency is 0.

#### Instructions
- Obtain the [yeast chromosome 1 sequence](http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr01.fsa)
If unable to access the file, get it from [here](https://github.com/ajodeh-juma/bixcop-2021-python/blob/main/modules/python/data/test/chr01.fsa)
- Write function(s) to compute the base frequencies.

#### Reporting
Report your output as:<br>
```A: 0.1```<br>```C: 0.2```<br>```G: 0.2```<br>```T: 0.5```
``````

