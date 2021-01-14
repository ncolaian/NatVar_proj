#! /usr/bin/env perl

# I will use this file to generate a metadata file contining all the genes within each fliC's genomic neighborhood (10kb)
# I will also use this to create individual fasta files for the neighborhood.


use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Path::Class;
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use Data::Dumper;

#My Variables
my $help = 0;
my $man = 0;
my $name_file; #Should be in the format of gene-genome
my $data_dir; #Directory of the dattabase
my $gff_extension = ".gff"; #Extension of the gene information
my $fasta_ext = ".faa";
my $out_dir;

#Read in the variables from the command line
GetOptions( 'man|m'   =>  \$man,
            'help|h'  =>  \$help,
            'flic_names|f=s' => \$name_file,
            'data_direct|dd=s' => \$data_dir,
            'gff_ext|g:s'  => \$gff_extension,
            'fasta_ext|e:s'   => \$fasta_ext,
            'out_dir|o=s'  => \$out_dir,
            ) || die("There was an error in the command line arguements\n");

# Pod Usage for the manual and the help pages
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose =>3) }

# Setup logging environment
my $logger = get_logger();

## Main ##

#Get the neighborhood files
open my $IN, "<", $name_file;
open my $OUT, ">", "$out_dir/defined_neighborhood.tsv";
mkdir "$out_dir/fasta_data";
#print $data_dir;
#exit 42;

while ( <$IN> ){
   chomp $_;
   my $gene = (split /-/, $_)[0];
   my $genome = (split /-/, $_)[1];
   
   my $neighborhood = get_start_and_end($gene, $genome);
   print $OUT "$_\t$gene\t$genome\t$neighborhood\n";
   
   make_fasta($gene, $genome, $neighborhood);
}

close $OUT;
close $IN;


## Subroutines ##

## This needs the gene name passed into it along with genome name
sub get_start_and_end {
   my ($name, $genome) = @_; 
   
   #I first will grep the line with the gene
   my $data_line = `grep $name $data_dir/$genome/$genome$gff_extension`;
   
   # I will then goo through and get the start and stop for the gene
   my @line_array = split /\t/, $data_line;
   my $gene_start = $line_array[3];
   my $gene_end = $line_array[4];
   
   #With this in mind, I will now grep the before and after 50 lines. This should narrow my search to make the algy faster
   my $before_lines = `grep -B 50 $name $data_dir/$genome/$genome$gff_extension`;
   my $after_lines = `grep -A 50 $name $data_dir/$genome/$genome$gff_extension`;
   
   my @before_lines_array = split /\n/, $before_lines;
   my @after_lines_array = split /\n/, $after_lines;
   
   #get rid of the original gene
   pop @before_lines_array;
   shift @after_lines_array;
   
   #holder for the names of the genes
   my $neighborhood_genes;
   
   #I now will go through the before lines and save the neighborhood genes
   my $keep = 0;
   foreach my $line ( @before_lines_array ) {
      if ( $keep == 1 ) {
         my $id = get_geneID_from_gff_line($line);
         $neighborhood_genes .= "$id,";
         next;
      }
      my $end = (split /\t/, $line)[4];
      #10kb neighborhood
      if ( $end >= ($gene_start - 10000) ) {
         $keep = 1;
         my $id = get_geneID_from_gff_line($line);
         $neighborhood_genes .= "$id,";
         next;
      }
      else {
         next;
      }
   }
   
   foreach my $line ( @after_lines_array ) {
      if ( $keep == 0 ) {
         next;
      }
      my $start = (split /\t/, $line)[3];
      #10kb neighborhood
      if ( $start >= ($gene_end + 10000) ) {
         $keep = 0;
         next;
      }
      else {
         my $id = get_geneID_from_gff_line($line);
         $neighborhood_genes .= "$id,";
      }
   }
   
   chop $neighborhood_genes; #This gets rid of the last ,
   
   return($neighborhood_genes);
}

sub get_geneID_from_gff_line {
   my ($line) = @_;
   
   my @liray = split /\t/, $line;
   my $id_portion = pop @liray;
   
   my @portions_ar = split /;/, $id_portion;
   my $id = (split /=/, $portions_ar[0])[1];
   return($id);
}

sub make_fasta {
   my ($gn, $gnm, $nb) = @_;
   open my $FASTA,">" ,"$out_dir/fasta_data/$gn-$gnm.faa";
   
   my @ngenes = split /,/, $nb;
   
   foreach my $gene ( @ngenes ) {
      my $data = `grep -A 1 $gene $data_dir/$gnm/$gnm$fasta_ext`; #This means that the fasta files need to be in a single line format
      print $FASTA $data;
   }
   close $FASTA;
}

__END__
=head1 TITLE



=head1 VERSION



=head1 INCLUDED MODULES

Getopt::Long;
Pod::Usage;
Carp;
Readonly;
Path::Class;
Data::Dumper;
Log::Log4perl qw(:easy);
Log::Log4perl::CommandLine qw(:all);

=head1 INHERIT

=head1 SYNOPSIS

=head1 PARAMETERS

=head1 CONFIGURATION AND ENVIRONMENT



=head1 DEPENDENCIES


    
=head1 INCOMPATIBILITIES

    None reported.

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests	
	
=head1 AUTHOR

Nicholas Colaianni
contact via C<< <ncolaian@live.unc.edu> >>

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2018, Nicholas Colaianni
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.

=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

=cut