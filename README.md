![JCAST](https://github.com/ed-lau/jcast/blob/master/images/logo2.png?raw=True)

# Junction Centric Alternative Splicing Translator v.0.3.1

JCAST (Junction Centric Alternative Splicing Translator) takes in alternative splicing events and returns custom protein sequence databases for isoform analysis.

## Getting Started

#### Requirements

Install Python 3.7+ and pip. See instructions on Python website for specific instructions for your operating system.

JCAST can be installed from PyPI via pip. We recommend using a virtual environment.

    $ pip install jcast

#### Running
	
Launch JCAST as a module (Usage/Help):
		
	$ python -m jcast

Alternatively:

    $ jcast

Example command: 
		
	$ python -m jcast  data/encode_human_pancreas/ data/gtf/Homo_sapiens.GRCh38.89.gtf data/gtf/Homo_sapiens.GRCh38.89.gtf data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa -o encode_human_pancreas
    
To test that the installation can load test data files in tests/data (sample rMATS file and human chr 15 genome files)

    $ pip install tox
    $ tox

To run JCAST using the test files and print the results to Desktop

    $ python -m jcast {j}/tests/data/rmats {j}/tests/data/genome/Homo_sapiens.GRCh38.89.chromosome.15.gtf  {j}/tests/data/genome/Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz -o ~/Desktop
    
where {j} is replaced by the path to JCAST.

## An example using JCAST to generate custom databases from ENCODE

The following is an example using ENCODE public RNA-seq dataset to generate a cardiac-specific database with JCAST.

#### Download RNA-Seq from ENCODE: 
As an example, we will download the .fastq files from ENCODE adult human heart [dataset 1](https://www.encodeproject.org/experiments/ENCSR436QDU/) and
     [dataset 2](https://www.encodeproject.org/experiments/ENCSR391VGU/).
     
#### Align the FASTQ files to a reference genome 
Read alignment can be done using STAR >= v.2.5.0, e.g.,:

	$ STAR --runThreadN 10 --genomeDir path/to/GRCh38/STARindex --sjdbGTFfile path/to/Homo_sapiens.gtf --sjdbOverhang 100 --readFilesIn ./ENCFF781VGS.fastq.gz ./ENCFF466ZAS.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./STAR_aligned/b1t1/
    $ STAR --runThreadN 10 --genomeDir path/to/GRCh38/STARindex --sjdbGTFfile path/to/Homo_sapiens.gtf --sjdbOverhang 100 --readFilesIn ./ENCFF731CDK.fastq.gz ./ENCFF429YOS.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./STAR_aligned/b2t1/

Note: Arguments including runThreadN and sjdbOverhang should be customized to suit your system and data files. Please refer to the STAR documentations for details.

#### Identify transcript splice junctions 
Splice junctions can be found using [rMATS](http://rnaseq-mats.sourceforge.net) with the .bam files following STAR. Please refer to the rMATS instructions for latest commands. The following
example was tested using rmats-turbo-0.1 running in Docker and using rMATS v.4.1.0/Python 3.7. Support for [stringtie](https://ccb.jhu.edu/software/stringtie/) assembled transcripts will be implemented in a future version.

#### Set up a Virtual Environment for rMATS turbo 0.1 in Python 2.7 (only if needed)

#### Install the rMATS image
Follow instructions from rMATS and docker specific to your OS. E.g.:

    $ sudo docker load -i rmats-turbo-0.1.tar

#### Prepare the /rMATS subdirectory 
Copy the individual .bam files from STAR into the rMATS subdirectory and rename them b1t1.bam, b1t2.bam, b2t1.bam, b2t2.bam, etc. Copy the GTF file from the Genomes folder as GRCm38.gtf. Write a b1.txt file with a text editor containing the following docker virtual directories:

    /data/b1t1.bam,/data/b1t2.bam
 
Write a b2.txt file

    /data/b2t1.bam,/data/b2t2.bam
 
Go back to the data directory and run the rMATS image. The -v flag mounts the host directory into the docker container at /data, which corresponds to the visual directories in the b1.txt and b2.txt files.

    $ sudo docker run -v path/to/data/directory:/data rmats:turbo01 --b1 /data/b1.txt --b2 /data/b2.txt --gtf /data/GRCh38.gtf --od /data/output -t paired  --nthread 4 --readLength 101 --anchorLength 1

Note: Arguments including nThread, readLength, and anchorLength should be customized to suit your system and data files. Please refer to the [rMATS](https://github.com/Xinglab/rmats-turbo) documentations for details.

Run the JCAST Python program specifying the path to the rMATS output directory, the genome sequence, as well as the GTF annotation file:
 
    $ python -m jcast path/to/rMATS/output/encode_human_heart/ path/to/gtf/Homo_sapiens.GRCh38.89.gtf path/to/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa -o encode_human_heart

#### FASTA output
JCAST outputs FASTA databases which can be further filtered and combined using any scripting languages, or can be used directly for database search
in virtually any shotgun proteomics database search engines (e.g., SEQUEST, Crux/Tide, Maxquant, MS-GF+).
    
JCAST may output the following FASTA files (note depending on the used settings and input files, not all FASTA files may be present):

* **xxx_canonical.fasta** -- This file contains protein sequences from splice junctions that are identical to SwissProt canonical sequences. The FASTA entries are named according to UniProt convention.    
* **xxx_T1.fasta** -- This file contains noncanonical sequences translated from splice junctions. Tier 1 junctions are translated in frame according to annotated GTF frames, did not encounter frameshift or premature stop codon, and are successfully joined back to full-length SwissProt sequences.
* **xxx_T2.fasta** -- Tier 2 junctions are translated in frame according to annotated GTF frames, did not encounter premature stop codon, and are successfully joined back to full-length SwissProt sequences, but have encountered a possible frameshift (length differences in exons not multiples of 3).
* **xxx_T3.fasta** -- Tier 3 junctions did not encounter premature stop codon, and are successfully joined back to full-length SwissProt sequences, but using a translation frame different from that annotated in the supplied GTF (they should be rare).
* **xxx_T4.fasta** -- Tier 4 junctions were forced-translated when one of the two alternative junction slices encountered a premature stop codon but could be translated using one of three frames into a peptide fragment at least a certain proportion in length as the successfully translated slice (see _params.py_). These sequences should be either excluded from database search or interpreted with a great amount of caution.
* **xxx_T#_orphan.fasta** -- These fragments were translated according to their tiers but could not be joined back to the canonical SwissProt sequence through the stitch length (see _params.py_ for defaults). These sequences should be either excluded from database search or interpreted with a great amount of caution.

Noncanonical FASTA entries have the following naming convention:

```
>sp|Q91VW5|GOGA4_MOUSE|ENSMUSG00000038708|MXE1|0|chr9|118560742:118560872|118565557:118565667|+2|r521|T1 sp|Q91VW5|GOGA4_MOUSE Golgin subfamily A member 4 OS=Mus musculus OX=10090 GN=Golga4 PE=1 SV=2
```

The vbar(|)-delimited parts denote the following:
1. Knowledgebase name, from canonical SwissProt protein entry (sp)
2. UniProt accession, from canonical SwissProt protein entry (Q91VW5)
3. UniProt name, from canonical SwissProt protein entry (GOGA4_MOUSE)
4. Annotated gene name (ENSMUSG00000038708)
5. rMATS junction type and order (MXE1)
6. Input file row name (0)
7. Chromosome (chr9)
8. Anchor exon start and end (118560742:118560872)
9. Alternative exon start and end (118565557:118565667)
10. Translated strand and phase (+2)
11. Minimal skipped junction count (sjc) in rMATS preceded by r (r521)
12. Tier (T1) 


### Dependencies

JCAST has been tested in Python 3.7, 3.8, 3.9 and uses the following packages:

```
biopython>=1.78
gtfparse>=1.2.1
pandas>=1.3.0
requests>=2.24.0
tqdm>=4.61.2
scikit-learn==0.24.2
matplotlib==3.4.2
```

## Additional Information

Additional details on troubleshooting and result interpretation can be found in our publication in [STAR Protocols](https://www.sciencedirect.com/science/article/pii/S2666166720301258).

## Contributing

Please contact us if you wish to contribute, and submit pull requests to us.


## Authors

* **Edward Lau, PhD** - *Code/design* - [ed-lau](https://github.com/ed-lau)
* **Maggie Lam, PhD** - *Code/design* - [Maggie-Lam](https://github.com/Maggie-Lam)


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
