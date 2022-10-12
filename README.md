# Junction Centric Alternative Splicing Translator v.0.3.5

JCAST (Junction Centric Alternative Splicing Translator) takes in alternative splicing events and returns custom protein sequence databases for isoform analysis.

See https://ed-lau.github.io/jcast/ for Documentation and Usage.

## Installation

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
		
	$ python -m jcast  data/encode_human_pancreas/ data/gtf/Homo_sapiens.GRCh38.89.gtf data/gtf/Homo_sapiens.GRCh38.89.gtf data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa -o encode_human_pancreas -q 0 1 -r 1 -m -c
    
To test that the installation can load test data files in tests/data (sample rMATS file and human chr 15 genome files)

    $ pip install tox
    $ tox

To run JCAST using the test files and print the results to Desktop

    $ python -m jcast {j}/tests/data/rmats {j}/tests/data/genome/Homo_sapiens.GRCh38.89.chromosome.15.gtf  {j}/tests/data/genome/Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz -o ~/Desktop


## All Arguments
```
python -m jcast -h
usage: __main__.py [-h] [-o OUT] [-r READ] [-m] [-c] [-q q_lo q_hi] [--g_or_ln G_OR_LN] rmats_folder gtf_file genome

jcast retrieves transcript splice junctionsand translates them into amino acid sequences

positional arguments:
  rmats_folder          path to folder storing rMATS output
  gtf_file              path to Ensembl gtf file
  genome                path to genome file

optional arguments:
  -h, --help            show this help message and exit
  -o OUT, --out OUT     name of the output files [default: psq_out]
  -r READ, --read READ  the lowest skipped junction read count for a junction to be translated [default: 1]
  -m, --model           models junction read count cutoff using a Gaussian mixture model [default: False]
  -c, --canonical       write out canonical protein sequence even if transcriptslices are untranslatable [default: False]
  -q q_lo q_hi, --qvalue q_lo q_hi
                        take junctions with rMATS fdr within this threshold [default: 0 1]
  --g_or_ln G_OR_LN     Switch on distribution to use for low end of histogram, 0 for Gamma, anything else for LogNorm


```


## Dependencies

JCAST has been tested in Python 3.7, 3.8, 3.9 and uses the following packages:

```
biopython>=1.78
gtfparse>=1.2.1
pandas>=1.3.0
requests>=2.24.0
tqdm>=4.61.2
scikit-learn>=1.0
matplotlib==3.4.2
scipy>=1.7.0
```

## Known Issues
* rMATS output with rows containing `NA` as gene name can fail.
* Upstream analyses should be performed using an unmasked genome. Currently JCAST cannot handle masked nucleotides (`N`).

## Additional Information

Additional details on troubleshooting and result interpretation can be found in our publication in [STAR Protocols](https://www.sciencedirect.com/science/article/pii/S2666166720301258).

## Contributing

Please contact us if you wish to contribute, and submit pull requests to us.


## Authors

* **Edward Lau, PhD** - *Code/design* - [ed-lau](https://github.com/ed-lau)
* **Maggie Lam, PhD** - *Code/design* - [Maggie-Lam](https://github.com/Maggie-Lam)
* **Robert Wes Ludwig, BSc** - *Modeling* - [WesLudwig](https://github.com/WesLudwig)


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
