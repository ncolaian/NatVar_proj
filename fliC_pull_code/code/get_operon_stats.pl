#! /usr/bin/env perl

# This file is the second script from neighborhood genes that grabs the genes within 10kB of each other. The idea is to find out ow many fliC operons there are in each genome
# The file to analyze will be the neighborhood tsv file.

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
my $neighborhood_file;
my $outfile;

#Read in the variables from the command line
GetOptions( 'man|m'   =>  \$man,
            'help|h'  =>  \$help,
            'neighbor|n=s'  => \$neighborhood_file,
            'outfile|o=s'  => \$outfile,
            ) || die("There was an error in the command line arguements\n");

# Pod Usage for the manual and the help pages
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose =>3) }

# Setup logging environment
my $logger = get_logger();

## Main ##
open my $IN, "<", $neighborhood_file;
open my $OUT, ">", $outfile;

print $OUT "Genome\t" . "num_genes\t" . "num_operons\t" . "count_of_small_operon\t" . "count_of_large_operon\t" . "num_in_operon\n";
my %done_genomes;
while ( <$IN> ) {
   chomp $_;
   my @split_line = split /\t/, $_;
   
   #check to see if the genome has been done yet
   if ( $done_genomes{$split_line[2]} ) {
      next;
   }
   $done_genomes{$split_line[2]} = 1;
   
   #get the operon info
   my ($operon_href, $num_genes) = get_operon_info($split_line[2]);
   
   #I will now print out the info
   my $num_operons = keys %$operon_href;
   
   print $OUT $split_line[2] . "\t$num_genes\t$num_operons";
   
   my $max = 0;
   my $min = 10;
   my @mem_counts;
   foreach my $operon (keys %$operon_href) {
      my $operon_members = scalar(split /,/, $operon_href->{$operon});
      push @mem_counts, $operon_members;
      if ( $operon_members > $max ) {
         $max = $operon_members;
      }
      if ( $operon_members < $min ) {
         $min = $operon_members;
      }
   }
   print $OUT "\t$min\t$max\t" . join("\t", @mem_counts) . "\n"
}

close($IN);
close($OUT);


## Subroutines ##
sub get_operon_info {
   my ($genome) = @_;
   
   #I will first get all the info for the genome
   my $genome_data = `grep $genome $neighborhood_file`;
   my @genome_file_by_line = split /\n/, $genome_data;
   
   my %genome_data_by_gene;
   my @genes;
   
   foreach my $line ( @genome_file_by_line ) {
      my @line_array = split /\t/, $line;
      #put the genome up to
      $genome_data_by_gene{$line_array[1]} = $line_array[3];
      push @genes, $line_array[1];
   }
   
   my $gene_number = scalar(@genes);
   
   my @sorted_genes = sort { $a <=> $b } @genes;
   
   my %operons;
   my $operon_number=1;
   for (my $i = 0; $i < $gene_number; $i++){
      if ( $i == 0 ) {
         $operons{1} = $sorted_genes[$i];
         next;
      }
      my $search_gene = $sorted_genes[$i];
      my $prev_gene_info = $genome_data_by_gene{$sorted_genes[$i-1]};
      if ( $prev_gene_info =~ /$search_gene/ ) {
         $operons{$operon_number} .= ",$search_gene";
      }
      else {
         $operon_number++;
         $operons{$operon_number} = $search_gene;
      }
   }
   
   return(\%operons, $gene_number);
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