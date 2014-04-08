#!/usr/bin/perl

###############################################################################
#         ____  _   _____   ___   ____  _____                      __         #
#        / __ \/ | / /   | |__ \ / __ \/ ___/___  ____ ___________/ /_        #
#       / /_/ /  |/ / /| | __/ // / / /\__ \/ _ \/ __ `/ ___/ ___/ __ \       #
#      / _, _/ /|  / ___ |/ __// /_/ /___/ /  __/ /_/ / /  / /__/ / / /       #
#     /_/ |_/_/ |_/_/  |_/____/_____//____/\___/\__,_/_/   \___/_/ /_/        #
#                                                                             #
###############################################################################
#                                                                             #
#       RNA2DSearch is an open source bioinformatics pipeline for searching   #
#       for short RNA secondary structures similar to user defined reference  #
#       structures on a genomic scale.                                        #
#                                                                             #
#       Contact: Russell.Hamilton@bioch.ox.ac.uk                              #
#                http://www.rna2dsearch.com                                   #
#                Department of Biochemistry,                                  #
#                University of Oxford, OX1 3QU                                #
#       Copyright (C) 2008 Russell S. Hamilton                                #
#                                                                             #
###############################################################################
# GNU Licence Details:                                                        #
# This program is free software; you can redistribute it and/or modify it     #
# under the terms of the GNU General Public License as published by the Free  #
# Software Foundation; either version 2 of the License, or (at your option)   #
# any later version.                                                          #
#                                                                             #
# This program is distributed in the hope that it will be useful, but WITHOUT #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       #
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    #
# more details (See http://www.opensource.org/licenses/gpl-license.php).      #
#                                                                             #
# You should have received a copy of the GNU General Public License along     #
# with this program; if not, write to the Free Software Foundation, Inc.,     #
# 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA                       #
###############################################################################

use strict;

my (@Databases, @QuerySeqs, $T, $L, $X, @SS_Seqs, );

#$T = 25;  # Folding Temperature
#$L = 64;  # Max basepair distance appart

$T = $ARGV[1];  # Folding Temperature
$L = $ARGV[2];  # Max basepair distance appart

chomp($ARGV[0]);
chomp($ARGV[1]);
chomp($ARGV[2]);

Run_RNALfold($ARGV[0]);


#-------------------------------------------------------------------------
sub Run_RNALfold {

my($Database, $i, $head, $length, $Input, $Query, $Output, @In, @In2, $j, 
   $Middle, $cnt, $Matches, $AllSeq, @Seq, );

$Database = $_[0];

($head,$length,$Input) = Deal_With_Input($Database);

$Output  = $Input;
$Output =~ s/.no_nl//g;
$Output .= ".output_$T" . "_$L";

`RNALfold -L $L -T $T < $Input > $Output`;

}

#-------------------------------------------------------------------------
sub Deal_With_Input {

my($InputFile, $Header, $Length, $WholeInput, $NewInputFile, );

$InputFile = $_[0];

open(INPUT, "$InputFile") || die "Can't open $InputFile: $!\n";
while (<INPUT>) { $WholeInput .= $_; }
close INPUT || die "Can't close $InputFile: $!\n";

# Deal with header
$WholeInput =~ m/(>.*)/;
$Header = $1;
chomp($Header);

# Remove Newline characters, for RNALfold input
$WholeInput =~ s/^\>.*//g;
$WholeInput =~ s/\n//g;

$Length = length($WholeInput);

$NewInputFile = $InputFile . ".no_nl";

open (NEWINPUT, ">$NewInputFile") || die "Can't open $NewInputFile: $!\n";
print NEWINPUT "$WholeInput\n";
close NEWINPUT || die "Can't close $NewInputFile: $!\n";

return ($Header, $Length, $NewInputFile);
}

#-------------------------------------------------------------------------
# FIN
#-------------------------------------------------------------------------
