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

my ($Databases, @QuerySeqs, $T, $L, @SS_Seqs, $Q,);

@QuerySeqs = 
   ("aaguaauuuucgugcucucaacaauugucgccgucacagauuguuguucgagccgaaucuuacu",     # GLS
    "ugcacaccucccucgucacucuugauuuuucaagagccuucgaucgaguaggugugca",           # ILS
    "ugcuggcaauccucuccaaaauacucgaaagaguauuucugcggagaguguugccaguac",         # G2LS
    "gugauacgcucgcggauuaacaaggggcucugaggagcuucggauauuccucagcgaucac"         # JLS
   );

@SS_Seqs   = 
   (".(((((..((((.((((.((((((((((.......))).)))))))...))))))))..)))))",     # GLS
    ".((((((((..((((((.(((((((....))))))).....)).)))).)))))))).",           # ILS
    "(((((((((..((((((.((((((((....))))))))....))))))..))))))))).",         # G2LS
    "(((((.((((.(.(((.......(((((((....)))))))......))).))))))))))"         #JLS
   );

# Command line arguments
$Databases = $ARGV[0];
$T         = $ARGV[1];  # Folding Temperature
$L         = $ARGV[2];  # Max basepair distance appart

Run_RNALfold($Databases, $Q);


#-------------------------------------------------------------------------
sub Run_RNALfold {

my($Database, $i, $head, $length, $Input, $Output, @Structures, $j, 
   $Middle, $cnt, $Matches, $AllSeq, @Seq, $OutputNL, @Dist_Score, 
   @Temp, $Dist_Out, $StructFile, $MatchesFile, $ContainsQ,    
   $ForesterFile, $raw, $hrs, $mins, $secs, $Stem, 
   $No_Q_BP, $No_R_BP,
   $Forest_Out_1, $Forest_Out_2, @Forest_Out_3, $Forest_Out_4, $Forest_Out_5,
   $LocalLength1, $LocalLength2,  
  );

$Database = $_[0];                  #Genomic seq being searched

print "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
print "+ Running $0\n";
print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

$Stem = int(rand(999));

# Read in the structures output file
print "+ Parsing structures file to an array... ";
open(OUTPUT, "$Database") || die "Can't open $Database: $!\n";
@Structures = <OUTPUT>;
close OUTPUT || die "Can't close $Database: $!\n";
print "DONE\n";
print "+    No. Structures = $#Structures\n";

$Database =~ s/.mfoldings.*//g;
$MatchesFile = "/usr/people/+user+/DB_MATCHES/" . $Database . ".mmatches_$T" . "_$L";
open(MATCHES,">$MatchesFile") || die "Can't open $MatchesFile: $!\n"; 

print MATCHES "+ DATABASE: $Database\n";
print MATCHES "+ $#Structures structures from RNALfold with L = $L\n\n";
print MATCHES "     1   5    10   15   20   25   30   35   40   45   50   55   60   65   70   75\n";
print MATCHES "     +...+....+....+....+....+....+....+....+....+....+....+....+....+....+....+\n\n";

print  "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

for($i=0; $i<=$#Structures; $i++)
   {
     #Deal With each structure output line from RNALfold
     $Structures[$i] =~ s/\( /\(/g;
     @Temp = split(/\s+/, $Structures[$i]);
     # ID
     $Temp[0] =~ s/ //g;       chomp($Temp[0]);
     # Structure
     $Temp[1] =~ s/ //g;       chomp($Temp[1]);
     # Energy
     $Temp[2] =~ s/\(//g; 
     $Temp[2] =~ s/\)//g;      chomp($Temp[2]);
     # Start
     $Temp[3] =~ s/ //g;       chomp($Temp[3]);
     # Sequence
     $Temp[4] =~ s/ //g;       chomp($Temp[4]);

     $Structures[$i] = "";  # hopefully free up a bit of memory

     for($j=0; $j<=$#QuerySeqs; $j++)
        {
          # Compare each structure to the query & get a score for it
          # Need to write struct & query to a file
          $StructFile = "/tmp/tempfile.MDS_$Stem\_$$\_$j";
          open(STRUCT, ">$StructFile") || die "Can't open $StructFile: $!\n";
          print STRUCT "$SS_Seqs[$j]\n$Temp[1]\n";
          close STRUCT || die "Can't close $StructFile: $!\n";

          $Dist_Out       =  `RNAdistance -DfhwcFHWC < $StructFile`;
          $Dist_Out       =~ m/f: (.*) h: (.*) w: (.*) c: (.*) F: (.*) H: (.*) W: (.*) C: (.*)/;
          $Dist_Score[$j] = "[$1]"; 
          $Dist_Score[$j] =~ s/ //g; 
          chomp($Dist_Score[$j]); 

          unlink($StructFile);
          undef($Dist_Out); 

          #RNAforester similarity
          $ForesterFile = "/tmp/tempfile.MDF_$Stem\_$$\_$j";

          $No_Q_BP = 2*($SS_Seqs[$j] =~ tr/(//); 
          $No_R_BP = 2*($Temp[1] =~ tr/(//);

          open (FOREST, ">$ForesterFile") || die "Can't open $ForesterFile: $!\n";
          print FOREST ">Query\n$QuerySeqs[$j]\n$SS_Seqs[$j]\n";
          print FOREST ">Match\n$Temp[4]\n$Temp[1]\n";
          close FOREST || die "Can't close $ForesterFile: $!\n";

          $Forest_Out_1 = `RNAforester --score --noscale -l < $ForesterFile`;
          chomp($Forest_Out_1); 
          $Forest_Out_1 =~ s/ //g;

          unlink($ForesterFile);
          unlink("rna.ps");
          unlink("x_1.ps");
          unlink("y_1.ps");
 
          $Dist_Score[$j] .= sprintf("[%.0f]", $Forest_Out_1);

          undef($Forest_Out_1); 
        }

     printf MATCHES "REP: ID=%s(%d) Length=%-3d Start=%-3d ",
                     $Temp[0], $cnt, length($Temp[1]), $Temp[3];
     printf MATCHES "Score=%s-%s-%s-%s Energy=[%s]\n",
                     $Dist_Score[0],$Dist_Score[1],$Dist_Score[2],$Dist_Score[3], $Temp[2];
     print  MATCHES "SEQ: $Temp[4]\n";
     print  MATCHES "SS_: $Temp[1]\n";
     $cnt++;
   }

undef($AllSeq);
undef(@Temp);

print MATCHES "\n";
print MATCHES "     +...+....+....+....+....+....+....+....+....+....+....+....+....+....+....+\n";
print MATCHES "     1   5    10   15   20   25   30   35   40   45   50   55   60   65   70   75\n\n";

close MATCHES || die "Can't close $MatchesFile: $!\n";

print  "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
print  "+ Summary:\n";
printf "+ Database entries              = %d\n", $#Structures;
printf "+ Total predicted structures    = %d\n", $#Structures+1; 
printf  "+ Total matching query in range = %d (%d %)\n", 
         $cnt, ($cnt/($#Structures+1)*100) ;

undef(@Structures);
undef($cnt);

print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
print "+ Finished Running $0\n"; 
print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
}

#----------------------------------------------------------------------------
sub MinMax {

my($answer);
                                                                                                       
if($_[2] =~ m/min/)
  {
    if($_[0] > $_[1]) { $answer = $_[1]; }
    else { $answer = $_[0]; }
  }
elsif($_[2] =~ m/max/)
  {
    if($_[0] > $_[1]){ $answer = $_[0]; }
    else { $answer = $_[1]; }
  }

return $answer;
}

#-------------------------------------------------------------------------
# FIN
#-------------------------------------------------------------------------
