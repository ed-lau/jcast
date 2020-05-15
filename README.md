# Junction Centric Alternative Splicing Translator v.0.1.0

JCAST (Junction Centric Alternative Splicing Translator) takes in alternative splicing events and returns custom protein sequence databases for isoform analysis.

## Getting Started

#### Requirements


Install Python 3.5+ and pip. See instructions on Python website for specific instructions for your operating system.

Set up virtual environment with venv at a designed path, e.g., ./venv/jcast
		
	$ python3.5 -m venv ./venv/jcast

Activate the virtual environment
	
		
	$ source ./venv/jcast/bin/activate

Check to ensure that the specific python build inside the venv is used
		
	
	$ which python

Install the packages in the requirements.txt file
		
	
	$ pip install -r requirements.txt


#### Running
	
Launch JCAST (Usage/Help)
		
	
	$ python3 jcast.py --help

Example command: 
		
	
	$ python jcast.py human data/encode_human_pancreas/ data/gtf/Homo_sapiens.GRCh38.89.gtf -o encode_human_pancreas
    
Deactivate the Virtual Environment upon completion
		
	
	$ deactivate


## Example using JCAST to generate custom databases from ENCODE

The following is an example using ENCODE public RNA-seq dataset to generate a cardiac-specific database with JCAST.

#### Download RNA-Seq from ENCODE: 
As an example, we will download the .fastq files from ENCODE adult human heart [dataset 1](https://www.encodeproject.org/experiments/ENCSR436QDU/) and
     [dataset 2](https://www.encodeproject.org/experiments/ENCSR391VGU/).
     
#### Align the FASTQ files to a reference genome 
Read alignment can be done using STAR v.2.5.0, e.g.,:


	$ STAR --runThreadN 10 --genomeDir path/to/STARindex --sjdbGTFfile path/to/Homo_sapiens.gtf --sjdbOverhang 100 --readFilesIn ./ENCFF781VGS.fastq.gz ./ENCFF466ZAS.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./STAR_aligned/b1t1/
    $ STAR --runThreadN 10 --genomeDir path/to/GRCh38/STARindex --sjdbGTFfile path/to/Homo_sapiens.gtf --sjdbOverhang 100 --readFilesIn ./ENCFF731CDK.fastq.gz ./ENCFF429YOS.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./STAR_aligned/b2t1/

#### Identify transcript splice junctions 
Splice junctions can be found using [rMATS](http://rnaseq-mats.sourceforge.net) with the .bam files following STAR. Please refer to the rMATS instructions for latest commands. The following
example was tested using rmats-turbo-0.1. Support for [stringtie](https://ccb.jhu.edu/software/stringtie/) assembled transcripts will be implemented in a future version.

#### Set up a Virtual Environment for rMATS in Python 2.7 (only if needed)

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

Run the JCAST Python program specifying the directory of the rMATS output as well as the GTF annotation file:

 
    $ python jcast.py human path/to/rMATS/output/encode_human_heart/ path/to/gtf/Homo_sapiens.GRCh38.89.gtf -o encode_human_heart

#### FASTA output
JCAST outputs FASTA databases which can be further filtered and combined using any scripting languages, or can be used directly for database search
in virtually any shotgun proteomics database search engines (e.g., SEQUEST, Crux/Tide, MS-GF+)
    
### Prerequisites

JCAST requires the following:

```
Python 3.5+
biopython
pandas
numpy
```


## Contributing

Please contact us if you wish to contribute, and submit pull requests to us.


## Authors

* **Edward Lau, PhD** - *Code/design* - [ed-lau](https://github.com/ed-lau)
* **Maggie Lam, PhD** - *Code/design* - [Maggie-Lam](https://github.com/Maggie-Lam)


## License

This project is licensed under the MIT License - see the [license.md](license.md) file for details


## Acknowledgments

* [PurpleBooth](https://github.com/PurpleBooth) for Github Readme template.