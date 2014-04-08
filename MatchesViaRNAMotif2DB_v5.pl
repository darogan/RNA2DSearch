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
use DBI;

my($DescFile, @Files2Run, $i, $Arm, $DB, $Temp, $Window, );

@Files2Run = $ARGV[0];
$DescFile  = $ARGV[1];   

for($i=0; $i<=$#Files2Run; $i++)
   {
      $Files2Run[$i] =~ m/dmel-(.*)-(.*)-r5.1.fasta.mmatches_(.*)_(.*)/;

      $Arm    = $1;
      $DB     = $2;
      $Temp   = $3;
      $Window = $4;

      print "+ Running: $Files2Run[$i]\n";
      print "  [Arm:$Arm DB:$DB, T:$Temp, Win:$Window]\n";

      CollateHits($Files2Run[$i]);

      print "+ Done: $Files2Run[$i]\n";
        
      undef($Arm); undef($DB); undef($Temp); undef($Window); 
   }

#----------------------------------------------------------------------------------------
sub CollateHits {

my(@Matches, $sql, $i, $j, @Results, @Hits, $num, @Details, $Annotation, $XSOME, $DBASE, 
   $TmpFile, $File, $SQL3, $SQL2, $CheckExists, $RNAMotif, $SQL, $TableDBName, 
   $dbh, $q, $row, $QuotedAnnotation, @SplitResults, @ResultsScore, @HighestScore, 
   $HighestScore, %Details, %Score, );

@Matches = ReadMatchesFile($_[0]);
print "  [Matches to be processed = $#Matches]\n";
$num=0;

$File = $_[0];
$File =~ m/^dmel-(.*)-(.*)-r5.1.fasta.mmatches*/;
$XSOME = $1;
$DBASE = $2;

for($i=0; $i<=$#Matches; $i++)
   {
#      $dbh = DBI->connect('DBI:mysql:RNA2DSearch_r5v1:192.168.4.127:3306',
#                          '+username+', '+password+',
#                          { RaiseError => 1, AutoCommit => 0} );

      $TmpFile = "/tmp/tmp.$DBASE.$XSOME.$Window.$$.$i";

      open(MOTIF,">$TmpFile.fasta") ;
      print MOTIF ">tmp\n$Matches[$i][1]\n";
      close MOTIF;

      @Results = `(rnamotif -descr $DescFile $TmpFile.fasta) 2>$TmpFile.error | rmprune | grep -v "^>" | grep -v "^#"`;

      if($#Results >= 0)
        {
          print "+\tFound Hit ($i,$num) ";

          $Details{'ID'} = {};
          $Matches[$i][0] =~ m/ID=(.*)\s+Length=(.*)\s+Start=(.*)\s+Score=(.*)\s+Energy=(.*)/;
          $Details{'ID'}     = $1; #ID
          $Details{'Length'} = $2; #Length
          $Details{'Start'}  = $3; #Start
          $Details{'Score'}  = $4; #Score
          $Details{'Energy'} = $5; #Energy
          $Details{'ID'} =~ s/\(.*//g;

	  $Details{'Score'} =~ m/\[(.*)\]\[(.*)\]-\[(.*)\]\[(.*)\]-\[(.*)\]\[(.*)\]-\[(.*)\]\[(.*)\]/;
          $Score{'GLS'}{'RNAdistance'}  = $1;
          $Score{'GLS'}{'RNAforester'}  = $2;
          $Score{'ILS'}{'RNAdistance'}  = $3;
          $Score{'ILS'}{'RNAforester'}  = $4;
          $Score{'G2LS'}{'RNAdistance'} = $5;
          $Score{'G2LS'}{'RNAforester'} = $6;
          $Score{'JLS'}{'RNAdistance'}  = $7;
          $Score{'JLS'}{'RNAforester'}  = $8;

          print "$Details{'ID'} ";

          $Annotation = GetGeneName($Details{'ID'});

          $RNAMotif = "";
          for($j=0; $j<=$#Results; $j++)
             {
                @SplitResults = split(/\s+/,$Results[$j]);
                $ResultsScore[$j] = sprintf("%.0f",$SplitResults[1]);


                chomp($Results[$j]);
                $RNAMotif .= "$Results[$j]xxx";
             }
          @HighestScore = sort {$a <=> $b} @ResultsScore;
          $HighestScore = $HighestScore[0];#$#HighestScore];

          print " Highest=$HighestScore ";
          $num++;
          $RNAMotif =~ s/xxx$//g;

          #DB Stuff
          $dbh = DBI->connect('DBI:mysql:RNA2DSearch_r5v1:192.168.4.127:3306',
                              '+username+', '+password+',
                              { RaiseError => 1, AutoCommit => 0} ); 


          $TableDBName = $DBASE;
          $TableDBName =~ s/-/_/g;

          #check exists
          $CheckExists = "";
          $sql = "SELECT Sequence, CG FROM RNA2DSearch_r5v1_$TableDBName " . 
                 "WHERE Sequence='$Matches[$i][1]' and CG='$Details{'ID'}'";
          $q   = $dbh->prepare($sql);
          $q->execute;
          while( $row = $q->fetchrow_hashref() )
               {
                 $CheckExists = "$row->{'Sequence'}_$row->{'CG'}";
               }     
 
          #update/select query Details[0] = GeneGC
          if( length($CheckExists) > 0 )
            {
              print "SEEN\n";
            }
          else
            {
              $QuotedAnnotation = $dbh->quote($Annotation);
              $SQL = "INSERT INTO RNA2DSearch_r5v1_$TableDBName VALUES(" . 
                     "'','$Matches[$i][1]','$Details{'ID'}','$XSOME','$Annotation'," .
                     "'$ARGV[0]'," . 
                     "'$Score{'GLS'}{'RNAdistance'}' ,'$Score{'GLS'}{'RNAforester'}' ," .
                     "'$Score{'ILS'}{'RNAdistance'}' ,'$Score{'ILS'}{'RNAforester'}' ," .
                     "'$Score{'G2LS'}{'RNAdistance'}','$Score{'G2LS'}{'RNAforester'}'," .
                     "'$Score{'JLS'}{'RNAdistance'}' ,'$Score{'JLS'}{'RNAforester'}' ," .
                     "'$ARGV[1]'," .
                     "'$HighestScore','$RNAMotif','$RNAMotif'," . 
                     "'','','','','','','',''," .
                     "0,'none','-') " ;
              $q = $dbh->prepare($SQL);
              $q->execute;

              print "NEW\n";
            }

          $dbh->commit;
          $dbh->disconnect;

          $Matches[$i] = "";

          undef($CheckExists);
          undef($SQL);
          undef($HighestScore);
          undef(@HighestScore);
          undef(@ResultsScore);
          undef($RNAMotif);
          undef(@Details);
          undef($Annotation);
          undef($QuotedAnnotation);
        } 
      
       unlink("$TmpFile.error"); 
       undef(@Results);
       unlink("$TmpFile.fasta");
#       $dbh->disconnect;
   }
undef(@Matches);
}

#----------------------------------------------------------------------------------------
sub ReadMatchesFile {
                                                                                                  
my ($MatchesFile, @MatchFile, $i, @Matches, $num, $num2, $num3, );
          
$MatchesFile = $_[0];                                                                                                  
open (MATCHES, "$MatchesFile") || die "Can't open $MatchesFile :$!\n";
@MatchFile = <MATCHES>;
close MATCHES || die "Can't close $MatchesFile :$!\n";
for($i=0; $i<=$#MatchFile; $i++)
   {
      chomp($MatchFile[$i]);
      if($MatchFile[$i] =~ m/^(REP: .*)/)
        {
          $Matches[$num][0] = $1;
          $Matches[$num][0] =~ s/REP: //g;
          $num++;
        }
      elsif($MatchFile[$i] =~ m/^SEQ: (.*)/)
        {
          $Matches[$num2][1] = $1;
          $Matches[$num2][1] = uc($Matches[$num2][1]);
          $Matches[$num2][1] =~ s/T/U/g;
          $num2++;
        }
      elsif($MatchFile[$i] =~ m/^(SS_: .*)/)
        {
          $Matches[$num3][2] = $1;
          $Matches[$num3][2] =~ s/SS_: //g;
          $num3++;
        }
   }                                                                                                  
undef(@MatchFile);
undef($MatchesFile);
undef($num);
undef($num2);
undef($num3);

return @Matches;
}

#------------------------------------------------------------------------------
sub GetGeneName {
                                                                                           
my ( $GeneLink, @Annot, $AnnotFile, $Annot, $Annotation, );
    
                                                                                           
$GeneLink = $_[0];

$AnnotFile = "/usr/people/+user+/DB_ANNOTATIONS/" .
             "fb_synonym_fb_2006_01_18.tsv";
                                                                                           
$Annot = `grep -i "$GeneLink" $AnnotFile`;
@Annot = split(/\s+/,$Annot);                   
                                                                        
if( $#Annot < 1)
  {
    $Annotation = "ANNOTATION MISSING FROM FILE\n";
  }
else
  {
    $Annotation = $Annot[1];
  }

undef($Annot);
                                                                                           
return $Annotation;
}

#----------------------------------------------------------------------------------------
# FIN
#----------------------------------------------------------------------------------------
