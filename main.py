"""
PSQ Junction Mapper

Usage:
    main.py <species> <rmats_folder> <gtf_file> <fasta_name>

Options:
    -h --help       Show this screen.
    -v --version    Show version.

Example:
    main.py mouse ~/rmats ~/Mus_musculus.GRCm38.85.gtf psq


"""

from get_jxn import Junction, Annotation
from get_seq import Sequence
from get_rma import RmatsResults


#
# Define output file name, overwrite existing file
#

def main(args):

    rmats_folder = args['<rmats_folder>']
    gtf_loc = args['<gtf_file>']
    out_file = args['<fasta_name>']
    species = args['<species>']

    #
    # Open the rMATS output file (MXE) here, rename the columns
    #

    rmats_results = RmatsResults(rmats_folder)

    #
    # Read the gtf file using the gtfpase package.
    # Then write as a pandas data frame.
    #
    gtf = Annotation(gtf_loc)
    gtf.read_gtf()


    #
    # Main loop through every line of each of the five rMATS files to make junction object, then translate them
    #

    for rmats_result in rmats_results.__dict__:

        rma = rmats_results.__getattribute__(rmats_result)
        for i in range(len(rma)):

            junction = Junction(id=rma.id[i],
                                gene_id=rma.gene_id[i],
                                strand=rma.strand[i],
                                gene_symbol=rma.gene_symbol[i],
                                chr=rma.chr[i],
                                anc_es=rma.anc_es[i],
                                anc_ee=rma.anc_ee[i],
                                alt1_es=rma.alt1_es[i],
                                alt1_ee=rma.alt1_ee[i],
                                alt2_es=rma.alt2_es[i],
                                alt2_ee=rma.alt2_ee[i],
                                down_es=rma.down_es[i],
                                down_ee=rma.down_ee[i],
                                junction_type=rma.jxn_type[i],
                                species=species,)

            print(str(i) + ' of ' + str(len(rma)))


            #
            # Subset the gtf file by the current gene_id
            #

            junction.get_translated_region(gtf)

            #
            # Get translation phase from GTF file.
            # If there is no phase found in the GTF, use phase -1 for now.
            # To do: look more closely into GTF file, or try translating from all frames
            #
            junction.get_translated_phase(gtf)

            #
            # Trim slice coordinates by translation starts and ends
            #
            junction.trim()

            #
            # Initiate a sequence object that copies most of the junction information
            #
            sequence = Sequence(junction)

            #
            # Get nucleotide sequences of all slices by REST API
            # (anchor, alternative-1, alternative-2, downstream)
            # Conjoin alternative exons to make slice 1 and 2,
            #
            sequence.make_slice()

            #
            # Translate into peptides
            #
            sequence.translate()

            #
            # Write into fasta
            #
            sequence.write_to_fasta(out_file)



from docopt import docopt

if __name__ == "__main__":
   args = docopt(__doc__, version='PSQJunctionMapper 0.1')
   print(args)
   main(args)

# python3 main.py rmats out_file

