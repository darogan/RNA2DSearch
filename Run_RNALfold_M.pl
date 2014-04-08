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
srand;
my (@Databases, @Headers, @Sequences, );

chomp($ARGV[0]); #DB file
chomp($ARGV[1]); #temp
chomp($ARGV[2]); #window

Run_RNALfold($ARGV[0],$ARGV[1],$ARGV[2]);

#-------------------------------------------------------------------------
sub Run_RNALfold {

my( $Database, $i, $SEQOUTFile, $Rand, $BIGOUTFile, $STRUCTOUTFile, $Seq, 
    @STRUCT_OUT, @Temp, $j, $T, $L, $num, $num2,  );

$Database = $_[0];
$T        = $_[1];
$L        = $_[2];

Deal_With_Input($Database);

$Rand          = int(rand(9999));
$BIGOUTFile    = "/usr/people/+user+/DB_FOLDINGS/" . $Database . ".mfoldings_$T" . "_$L";
$SEQOUTFile    = "/tmp/tempfile.MSeq.$Rand.$$";
$STRUCTOUTFile = "/tmp/tempfile.MStr.$Rand.$$";

printf "No. Sequences = %d\n", $#Sequences+1;


open(BIGOUT, ">$BIGOUTFile") || die "Can't open $BIGOUTFile: $!\n";


$num2=0;
for($i=0; $i<=$#Headers; $i++)
   {
     $num=0; 

     # Print each seq to file for RNALfold
     open(SEQOUT, ">$SEQOUTFile") || die "Can't open $SEQOUTFile: $!\n";
     print SEQOUT "$Sequences[$i]\n";
     close SEQOUT || die "Can't close $SEQOUTFile: $!\n";    

     # Run RNALfold
     `RNALfold -L $L -T $T < $SEQOUTFile > $STRUCTOUTFile`;

     # Deal with individual structure predictions
     open(STRUCTOUT, "$STRUCTOUTFile") || die "Can't open $STRUCTOUTFile: $!\n";
     @STRUCT_OUT = <STRUCTOUT>;
     close STRUCTOUT || die "Can't close $STRUCTOUTFile: $!\n";

     for($j=0; $j<=$#STRUCT_OUT; $j++)
        {
          if($STRUCT_OUT[$j] =~ m/^[\.\(\)]/)
            {
              #Deal With each structure output line from RNALfold
              $STRUCT_OUT[$j] =~ s/\( /\(/g;
              @Temp = split(/\s+/, $STRUCT_OUT[$j]);
              # Structure
              $Temp[0] =~ s/ //g;       chomp($Temp[0]);
              # Energy
              $Temp[1] =~ s/\(//g;
              $Temp[1] =~ s/\)//g;      chomp($Temp[1]);
              # Start
              $Temp[2] =~ s/ //g;       chomp($Temp[2]);

              # Get the sequence for the structure
              $Seq = substr($Sequences[$i],$Temp[2]-1,length($Temp[0]));

              # Print to big output file
              print BIGOUT "$Headers[$i] $Temp[0] $Temp[1] $Temp[2] $Seq\n";
              $num++;
              $num2++;
            }
        }
     printf "Dealt with %4d has %5d structures from seq length %6d\n", 
            $i+1, $num, length($Sequences[$i]);
   }

close BIGOUT || die "Can't close $BIGOUTFile: $!\n";

unlink($SEQOUTFile);
unlink($STRUCTOUTFile);

print "Finished with $num2 structures in total\n";

}

#-------------------------------------------------------------------------
sub Deal_With_Input {

my($InputFile, @whole, $i, $header_count, $seq_count,  );

$InputFile = $_[0];

open(INPUT, "$InputFile") || die "Can't open $InputFile: $!\n";
@whole = <INPUT>;
close INPUT || die "Can't close $InputFile: $!\n";

for($i=0; $i<=$#whole; $i++)
   {
     $whole[$i] = uc($whole[$i]);

     if($whole[$i] =~ />/)  
       {
         $Headers[$header_count] = $whole[$i];
         chomp($Headers[$header_count]);
         $seq_count++;
         $header_count++;
       }
     else                      
       {
         $Sequences[($seq_count-1)] = $Sequences[$seq_count-1] . $whole[$i];
       }
   }

for($i=0; $i<=$#Sequences; $i++)
   {
     @Sequences[$i] =~ s/\n//g;
     @Headers[$i]   =~ s/ .*//g;
     @Headers[$i]   =~ s/>//g;
   }

undef(@whole);
undef($seq_count);
undef($header_count);
}

#-------------------------------------------------------------------------
# FIN
#-------------------------------------------------------------------------
