BEASTGen Tool

This package contains a simple program for transforming data from one format to another.

The program is run using the following command line:

beastgen <template_filename> <input_filename> <output_filename>

For example:

beastgen to_fasta.template dengue.nex dengue.fasta

This will change the NEXUS file into a FASTA file by using the to_fasta template (a FASTA
to NEXUS template is also given in the templates folder).

Also in the templates folder is 'beast_example.template' which demonstrates how to
set up a template from a BEAST XML file.

beastgen -chain_length 10000000 -log_every 1000 beast_example.template dengue.nex

More detailed instructions are to be found on the BEAST website:

http://beast.bio.ed.ac.uk/
