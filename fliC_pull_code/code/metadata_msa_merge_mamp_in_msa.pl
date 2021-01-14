#! /usr/bin/env perl

# This is a custom script to pull out the flagellin sequence from a msa and link to metadata

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
my $mamp_res;
my $metadata;
my $out_dir;
my $res_msa; #this is the msa created from getting the first 120 aa form flagellin genes

#Read in the variables from the command line
GetOptions( 'man'   =>  \$man,
            'help'  =>  \$help,
            'metadata|md=s'  =>  \$metadata,
            'out_dir|o=s'   =>  \$out_dir,
            'res_msa|r=s'   =>  \$res_msa,
            ) || die("There was an error in the command line arguements\n");

# Pod Usage for the manal and the help pages
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose =>3) }

# Setup logging environment
my $logger = get_logger();

## Main ##
open my $META, "<", $metadata || die "Can't open metadata file\n";
my $header = <$META>;

my %id_hash;
while (<$META>) {
    chomp $_;
    my @line = split /\t/, $_; 
    
    #need to get rid of random commas in the metadata
    $_ =~ s/,//g;
    
    #save the metadata line with the id
    $_ =~ s/\t/,/g;
    $id_hash{$line[0]} = $_;
}
close $META;

open my $KOUT, ">", "$out_dir/known_mamp_seq_metadata.csv";
open my $BOUT, ">", "$out_dir/bad_full_alignment.txt";
open my $MSA, "<", $res_msa;

chomp $header;
$header =~ s/\t/,/g;
print $KOUT $header . ",gene-genome";
#add header for full mamp seq
for ( my $i = 1; $i <= 22; $i++){
    print $KOUT ",pos$i";
}
print $KOUT "\n";

my $flg_seq = `grep -A 1 flg $res_msa`;
my @lines = split /\n/, $flg_seq;
my @seq_array = split //, (split /\n/, $flg_seq)[1];
my @positions;
for (my $i = 0; $i < scalar(@seq_array);$i++) {
    if ( $seq_array[$i] !~ /-/ ) {
        push @positions, $i;
    }
    else {
      if ( scalar(@positions > 0) ) {
         print "WARNING: The alignment of the MAMP is not continuous for the first 14 amino acids\n";
      }
    }
    if ( scalar(@positions == 14) ) {
      last;
    }
}

my $id;
while ( <$MSA> ) {
    chomp $_;
    if ( $_ =~ /^>/ ) {
        ($id = $_) =~ s/>//g; #doesnt change $_ but creates id equal to $_ w/o >
    }
    else {
        my @full_seq = split //, $_;
        my $seq;
        #start position from the alignment and then go up 22 spaces
        for ( my $i = $positions[0]; $i < scalar(@full_seq); $i++ ) {
            if ( $full_seq[$i] !~ /-/ ) {
               $seq .= $full_seq[$i];
            }
            if ( $seq && length($seq) == 22 ) {
               last;
            }
        }
        if ( (split //, $seq)[13] ne "D" ) {
            print "$id\'s sequence does not have a D at position 14. Consider looking at the sequence\n";
            print $BOUT "$id\n$seq\n";
            next;
        }
        
        if ( $seq =~ /X/ ) {
            print "$id\'s sequence has an X in the flg22 portion and will be removed.\n";
            print $BOUT "$id\n$seq\n";
            next;
        }
        my $genome = (split /-/, $id)[1];
        my @f_seq = split //, $seq;
        if ( $id_hash{$genome} ) {
            print $KOUT $id_hash{$genome} . "," . "$id," . join(",",@f_seq) . "\n";
        }
    }
}

close( $MSA );
close( $KOUT );
close( $BOUT );



## Subroutines ##

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

    Need to load muscle and BLAST

=head1 DEPENDENCIES

    muscle
    BLAST
    Master_aln
    
=head1 INCOMPATIBILITIES

    None reported.

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests	
	
=head1 AUTHOR

Nicholas Colaianni
contact via C<< <ncolaian@live.unc.edu> >>

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2017, Nicholas Colaianni
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