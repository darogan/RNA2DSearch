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
use Time::Local;

my ($Databases, $time_taken, @QuerySeqs, $T, $L, $X, @SS_Seqs, $SM, $Q,);

@QuerySeqs = 
   ("aagtaattttcgtgctctcaacaattgtcgccgtcacagattgttgttcgagccgaatcttact",
    "tgcacacctccctcgtcactcttgatttttcaagagccttcgatcgagtaggtgtgca");
@SS_Seqs   = 
   (".(((((..((((.((((.((((((((((.......))).)))))))...))))))))..)))))", 
    ".((((((((..((((((.(((((((....))))))).....)).)))).)))))))).");

# Command line arguments
$Databases = $ARGV[0];
$T         = $ARGV[1];  # Folding Temperature
$L         = $ARGV[2];  # Max basepair distance appart
$X         = $ARGV[3];  # Score cut off for results (grk_vs_ifactor = 24)
$SM        = $ARGV[4];  # scoring method for RNAdistance (f h w c F H C W)
$Q         = $ARGV[5];  # query id

$time_taken = time;
Run_RNALfold($Databases, $Q);


#-------------------------------------------------------------------------
sub Run_RNALfold {

my($Database, $i, $head, $length, $Input, $Query, $Output, @Structures, $j, 
   $Middle, $cnt, $Matches, $AllSeq, @Seq, $OutputNL, $Dist_Score, 
   @Temp, $Dist_Out, $StructFile, $MatchesFile, $ContainsQ, %ScoreFrequency,   
   $raw, $hrs, $mins, $secs, $Stem, 
  );

$Database = $_[0];                  #Genomic seq being searched
$Query    = $_[1];                  #Query seq & ss id for arrays above
$Output   = $Database . ".output_$T" . "_$L";  #Output file for structures

print "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
print "+ Running $0\n";
print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

$Stem = $Database;
$Stem =~ s/_.*//g;

#Check to see if the database file has already been structure predected
if( (FileCheck("$Output") == 0) || (FileCheck("$Database.no_nl") == 0) )
  {
    `perl Run_RNALfold.pl $Database $T $L $X`;
    print "+  Running RNALfold    ... ";
  }
elsif( (FileCheck("$Output") == 1) && (FileCheck("$Database.no_nl") == 1) )
  {
    print "+ Already Run RNALfold...";
  }
print " DONE\n";

# Read in the structures output file
print "+ Parsing structures file to an array... ";
open(OUTPUT, "$Output") || die "Can't open $Output: $!\n";
@Structures = <OUTPUT>;
close OUTPUT || die "Can't close $Output: $!\n";
print "DONE\n";
print "+    No. Structures = $#Structures\n";

# read in the no new line version of the sequence
print "+ Parsing database file to remove nl characters... ";
$OutputNL = $Database . ".no_nl";
open(OUTPUTNL, "$OutputNL") || die "Can't open $OutputNL: $!\n";
$AllSeq = <OUTPUTNL>;
close OUTPUTNL || die "Can't close $OutputNL: $!\n";
print "DONE\n";

$length = length($AllSeq);

print "+    No. nucleotides = $length\n";

$MatchesFile = $Database . ".matches_$T" . "_$L" . "_$X" . "_$SM" . "_$Q";;
open(MATCHES,">$MatchesFile") || die "Can't open $MatchesFile: $!\n"; 

print MATCHES "+ DATABASE: $Database\n";
print MATCHES "+ $#Structures structures from RNALfold with L = $L and $length nt\n\n";
print MATCHES "     1   5    10   15   20   25   30   35   40   45   50   55   60   65   70   75\n";
print MATCHES "     +...+....+....+....+....+....+....+....+....+....+....+....+....+....+....+\n\n";

print  "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

$| = 1;  # This disables buffering on stdout
my @Symbol = ("|", "/", "-", "\\", "|", "/", "-", "\\");
my $num=0;

#Draw a schematic for positions of matching sequences on database seq

for($i=0; $i<=$#Structures; $i++)
   {
     if($Structures[$i] =~ m/^[\.\(\)]/)
       {
         #Deal With each structure output line from RNALfold
         $Structures[$i] =~ s/\( /\(/g;
         @Temp = split(/\s+/, $Structures[$i]);
         # Structure
         $Temp[0] =~ s/ //g;       chomp($Temp[0]);
         # Energy
         $Temp[1] =~ s/\(//g; 
         $Temp[1] =~ s/\)//g;      chomp($Temp[1]);
         # Start
         $Temp[2] =~ s/ //g;       chomp($Temp[2]);

         # Compare each structure to the query & get a score for it
         # Need to write struct & query to a file
         $StructFile = "tempfile.A_$Stem";
         open(STRUCT, ">$StructFile") || die "Can't open $StructFile: $!\n";
         print STRUCT "$SS_Seqs[$Query]\n$Temp[0]\n";
         close STRUCT || die "Can't close $StructFile: $!\n";

         $Dist_Out   =  `RNAdistance -D$SM < $StructFile`;
         $Dist_Out   =~ m/$SM: (.*)/;
         $Dist_Score =  $1; chomp($Dist_Score); $Dist_Score =~ s/ //g;

         $ScoreFrequency{$Dist_Score}++;

         if($Dist_Score <= $X )
           {
             # Print outputs = overlaying progress report
#             printf "\r+ Progress[%s]: %d (%3.0f pcnt) %d (%3.0f pcnt)", 
#                     $Symbol[$num], 
#                     $cnt+1, 
#                     (($cnt+1) / $i) * 100,  
#                     $i+1, 
#                     ((($i+1)/$#Structures)*100+1);

             if($num == $#Symbol) { $num = 0; }
             else{ $num++; }

             # Find the subsequence associated with the structure
             $Seq[$i] = substr($AllSeq,$Temp[2]-1,length($Temp[0]));

             $Middle = substr($QuerySeqs[$Query], 
                              sprintf("%.0f",(length($QuerySeqs[$Query])/8)*3),
                              sprintf("%.0f",(length($QuerySeqs[$Query])/4)) );

             if($Seq[$i] =~ m/(.*)($Middle)(.*)/)
               {
                 $ContainsQ = "Y";
               }
             else
               {
                 $ContainsQ = "N";
               }


             printf MATCHES "REP: ID=%d(%d) Length=%-3d Start=%-3d ",
                             $i+1, $cnt, length($Temp[0]), $Temp[2];
             printf MATCHES "Score=%-3d Energy=[%s] Query=%s\n",
                             $Dist_Score, $Temp[1], $ContainsQ;
             print  MATCHES "SEQ: $Seq[$i]\n";
             print  MATCHES "SS_: $Temp[0]\n";
             $cnt++;
           }
       }
   }

undef($AllSeq);
undef(@Temp);

print MATCHES "\n";
print MATCHES "     +...+....+....+....+....+....+....+....+....+....+....+....+....+....+....+\n";
print MATCHES "     1   5    10   15   20   25   30   35   40   45   50   55   60   65   70   75\n\n";

close MATCHES || die "Can't close $MatchesFile: $!\n";

print  "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
print  "+ Summary:\n";
printf "+ Database sequence length      = %d\n", $length;
printf "+ Maximum possible structures   = %d\n", $length - $L + 1; 
printf "+ Total predicted structures    = %d (%d %)\n", 
        $#Structures+1, (($#Structures+1)/$length)*100; 
printf  "+ Total matching query in range = %d (%d %)\n", 
         $cnt, ($cnt/($#Structures+1)*100) ;

undef(@Structures);
undef($length);
undef($cnt);

#my $StatsFile = "$Database.stats_$T" . "_$L" . "_$X";
#open(STATS,">$StatsFile") || die "Can't open $StatsFile: $!\n";
#my $key;
#foreach $key (keys %ScoreFrequency)
#  {
#    print STATS "$key\t$ScoreFrequency{$key}\n";
#  }
#close STATS || die "Can't close $StatsFile: $!\n";
#
#undef(%ScoreFrequency);
#undef($key);
#
#`sort -n $StatsFile > $StatsFile.sorted_$T-$L-$X`;
#`Gnuplot.shell $StatsFile.sorted_$T-$L-$X &`;
#
#print "+ Stats file: $StatsFile.sorted_$T-$L-$X\n";

$time_taken = time - $time_taken;
$raw        = $time_taken;

$secs       = $time_taken % 60;
if(length($secs) == 1){$secs = "0" . $secs;}
$time_taken = ($time_taken - $secs) / 60;
$mins       = $time_taken % 60;
if(length($mins) == 1){$mins = "0" . $mins;}
$time_taken = ($time_taken - $mins) / 60;
$hrs        =  $time_taken % 24;
if(length($hrs) == 1){$hrs = "0" . $hrs;}

print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
print "+ Time = $hrs hrs: $mins mins: $secs secs (Raw=$raw)\n";
print "+ Finished Running $0\n"; 
print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
}

#----------------------------------------------------------------------------
sub FileCheck {

my($i, $FileToCheck, @DirFiles, $Presence,  );

$FileToCheck = $_[0];
@DirFiles    = <*>;

$Presence = 0;

for($i=0; $i<=$#DirFiles; $i++)
   {
     if($DirFiles[$i] =~ /\b$FileToCheck\b/)
       {
         $Presence = 1;
         last;
       }
   }

if($Presence == 1)
  {
    return 1;
  }
else
  {
    return 0;
  }
}

#-------------------------------------------------------------------------
# FIN
#-------------------------------------------------------------------------
