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
                                                                                           
system("clear");

my(%HoH, $number, @Unique, $i, $j, @Headers, @Sequences, @Dir, @XSome, @File, );

$Dir[0] = "/usr/people/+user+/Databases/DB_PROTEINCODING";
$Dir[1] = "/usr/people/+user+/Databases/DB_TRANSCRIPT";
$Dir[2] = "/usr/people/+user+/Databases/DB_CDS";

$File[0] = "_protein-coding-gene_dmel_RELEASE3-1.FASTA";
$File[1] = "_transcript_dmel_RELEASE3-1.FASTA";
$File[2] = "_CDS_dmel_RELEASE3-1.FASTA";

$XSome[0] = "2L";
$XSome[1] = "2R";
$XSome[2] = "3L";
$XSome[3] = "3R";
$XSome[4] = "X";

my $file_to_run;

for($i=0; $i<=$#Dir; $i++)
   {
     for($j=0; $j<=$#XSome; $j++)
        {
          $file_to_run = $Dir[$i] . "/" . $XSome[$j] . $File[$i];
          print "START : $file_to_run\n";
          Deal_With_Input($file_to_run);
#          print "+-+$#Sequences - $#Headers+-+\n";
          Parse_Out($i);

print "FIN   : $file_to_run\n";
#          print "$#Sequences - $#Headers $file_to_run\n";
          undef($file_to_run);
          undef(@Headers);
          undef(@Sequences);
        }
   }

print "\n\nTotal = $number\n";

TraverseHash();

#------------------------------------------------------------------------------
sub UniqueList {

my($CG, $ID, $FB, $DB, $SE, );

$CG = $_[0];
$ID = $_[1];
$FB = $_[2];
$DB = $_[3];
$SE = $_[4];

# keep longest sequences!

if($HoH{$CG})
  {
    print "+ ";
  }
else
  {
    print "- ";
  }

$HoH{$CG}{"FB_ID"} = $FB;
$HoH{$CG}{"GeneName"} = $ID;

if($DB == 0)
  {
    if( (length($SE >= $HoH{$CG}{"PC_Seq"})) )
      {
        $HoH{$CG}{"PC_Seq"} = $SE;
      }
  }
elsif($DB == 1)
  {
    if( (length($SE >= $HoH{$CG}{"TS_Seq"})) )
      {
         $HoH{$CG}{"TS_Seq"} = $SE;
      }
  }
elsif($DB == 2)
  {
    if( (length($SE >= $HoH{$CG}{"CDS_Seq"})) )
      {
        $HoH{$CG}{"CDS_Seq"} = $SE;
      }
  }


printf "%-20s %-20s %-20s %d %d\n", $CG, $ID, $FB, $DB, length($SE);

$number++;
}

#------------------------------------------------------------------------------
sub TraverseHash {

my($CGEntry, $role, $cnt, $num, $dbh, $q, $sql, $QuotedGeneName,);

$dbh = DBI->connect('DBI:mysql:RNALfold_Results:192.168.4.127:3306', '+username+', '+password+',
                    { RaiseError => 1, AutoCommit => 1} );
$dbh->{LongReadLen} = 16384;
$dbh->{LongTruncOk} = 1;


$num = 0;
$cnt = 0;

foreach $CGEntry ( keys %HoH ) 
   {
     if( ( (length($HoH{$CGEntry}{'PC_Seq'})  > 0) ||
           (length($HoH{$CGEntry}{'TS_Seq'})  > 0) || 
           (length($HoH{$CGEntry}{'CDS_Seq'}) > 0) ) &&
           (length($HoH{$CGEntry}{'PC_Seq'}) >= length($HoH{$CGEntry}{'TS_Seq'})  ) &&
           (length($HoH{$CGEntry}{'TS_Seq'}) >= length($HoH{$CGEntry}{'CDS_Seq'}) )    
       )
       {
         printf "%-8s: ", $CGEntry;

         printf "FB_ID=%-11s, ", $HoH{$CGEntry}{'FB_ID'};
         printf "GeneName=%-20s, ", $HoH{$CGEntry}{'GeneName'};

         printf "PC=%-6d TS=%-6d CDS=%-6d\n", length($HoH{$CGEntry}{'PC_Seq'}),
                                              length($HoH{$CGEntry}{'TS_Seq'}),
                                              length($HoH{$CGEntry}{'CDS_Seq'});

         $QuotedGeneName = $dbh->quote($HoH{$CGEntry}{'GeneName'});

         $sql = "INSERT INTO GeneInfo VALUES ('$CGEntry', '$HoH{$CGEntry}{'FB_ID'}', $QuotedGeneName, '$HoH{$CGEntry}{'PC_Seq'}', '$HoH{$CGEntry}{'TS_Seq'}','$HoH{$CGEntry}{'CDS_Seq'}')";

         $q = $dbh->prepare($sql);
         $q->execute;
         $q->finish;
         undef($q);

         $num++;
       }
     else
       {
         print "DODGY! - ";
         printf "%-8s: ", $CGEntry; 
         printf "PC=%-6d TS=%-6d CDS=%-6d\n", length($HoH{$CGEntry}{'PC_Seq'}),
                                              length($HoH{$CGEntry}{'TS_Seq'}),
                                              length($HoH{$CGEntry}{'CDS_Seq'});
         $cnt++;
       }
   }

$dbh->disconnect;

print "There are $num unique CG numbers ($cnt DODGY!)\n";
}

#------------------------------------------------------------------------------
sub Parse_Out {

my($i, $num, $j, $k, $tmp, @CG, @FB, @ID, $File_In,);
$File_In = $_[0];

$num=0;

for($i=0; $i<=$#Headers; $i++)
   {
     if( ($Headers[$i] !~ m/FLAGGED_PROBLEMATIC/) && 
         ($Headers[$i] =~ m/gene symbol:/) &&
         ($Headers[$i] !~ m/synonyms:/) &&
         ($Headers[$i] !~ m/gene_properties:/) && 
         ($Headers[$i] !~ m/snRNA symbol:/) && 
         ($Headers[$i] !~ m/snoRNA symbol:/) &&
         ($Headers[$i] !~ m/pseudogene symbol:/) &&
         ($Headers[$i] !~ m/tRNA symbol:/))
       {
          $Headers[$i] =~ m/(CG[0-9]+) gene symbol:([^ ]+) (FB[a-zA-Z0-9]+) /;
          $CG[$num] = $1;
          $ID[$num] = $2;
          $FB[$num] = $3;
#          print "$i\t\t$CG[$num] = $ID[$num] = $FB[$num] syn- prop- rna- pseudo- tRNA- snoRNA-\n";
          $num++;

          UniqueList($1,$2,$3,$File_In,$Sequences[$i]);
       }
     elsif( ($Headers[$i] !~ m/FLAGGED_PROBLEMATIC/) &&
         ($Headers[$i] !~ m/synonyms:/) &&
         ($Headers[$i] !~ m/gene_properties:/) &&
         ($Headers[$i] =~ m/snoRNA symbol:/) &&
         ($Headers[$i] !~ m/pseudogene symbol:/) &&
         ($Headers[$i] !~ m/tRNA symbol:/))
       {
          $Headers[$i] =~ m/(C[A-Z][0-9]+) snoRNA symbol:([^ ]+) (FB[a-zA-Z0-9]+) /;
          $CG[$num] = $1;
          $ID[$num] = $2;
          $FB[$num] = $3;
#          print "$i\t\t$CG[$num] = $ID[$num] = $FB[$num] syn- prop- rna- pseudo- tRNA- snoRNA+\n";
          $num++;

          UniqueList($1,$2,$3,$File_In,$Sequences[$i]);
       }
     elsif( (($Headers[$i] =~ m/FLAGGED_PROBLEMATIC/) ||
         ($Headers[$i] =~ m/synonyms:/) ||
         ($Headers[$i] =~ m/gene_properties:/) )&&
         ($Headers[$i] =~ m/snoRNA symbol:/) &&
         ($Headers[$i] !~ m/pseudogene symbol:/) &&
         ($Headers[$i] !~ m/tRNA symbol:/))
       {
          $Headers[$i] =~ m/(C[A-Z][0-9]+) snoRNA symbol:([^ ]+) .* (FB[a-zA-Z0-9]+) /;
          $CG[$num] = $1;
          $ID[$num] = $2;
          $FB[$num] = $3;
#          print "$i\t\t$CG[$num] = $ID[$num] = $FB[$num] syn- prop- rna- pseudo- tRNA- snoRNA+\n";
          $num++;

         UniqueList($1,$2,$3,$File_In,$Sequences[$i]);
       }

     elsif( ($Headers[$i] !~ m/FLAGGED_PROBLEMATIC/) &&
         ($Headers[$i] !~ m/synonyms:/) &&
         ($Headers[$i] !~ m/gene_properties:/) &&
         ($Headers[$i] !~ m/snRNA symbol:/) &&
         ($Headers[$i] !~ m/pseudogene symbol:/) &&
         ($Headers[$i] =~ m/tRNA symbol:/))
       {
          $Headers[$i] =~ m/(C[A-Z][0-9]+) tRNA symbol:([^ ]+) (FB[a-zA-Z0-9]+) /;
          $CG[$num] = $1;
          $ID[$num] = $2;
          $FB[$num] = $3;
#          print "$i\t\t$CG[$num] = $ID[$num] = $FB[$num] syn- prop- rna- pseudo- tRNA+\n";
          $num++;

           UniqueList($1,$2,$3,$File_In,$Sequences[$i]);
       }
     elsif( ($Headers[$i] !~ m/FLAGGED_PROBLEMATIC/) &&
         (($Headers[$i] =~ m/synonyms:/) || ($Headers[$i] =~ m/gene_properties:/) ) &&
         ($Headers[$i] !~ m/snRNA symbol:/) &&
         ($Headers[$i] !~ m/pseudogene symbol:/) &&
         ($Headers[$i] =~ m/tRNA symbol:/))
       {
          $Headers[$i] =~ m/(C[A-Z][0-9]+) tRNA symbol:([^ ]+) .* (FB[a-zA-Z0-9]+) /;
          $CG[$num] = $1;
          $ID[$num] = $2;
          $FB[$num] = $3;
#          print "$i\t\t$CG[$num] = $ID[$num] = $FB[$num] syn- prop+ rna- pseudo- tRNA+ \n";
          $num++;

         UniqueList($1,$2,$3,$File_In,$Sequences[$i]);
       }
     elsif( ($Headers[$i] !~ m/FLAGGED_PROBLEMATIC/) &&
         ($Headers[$i] !~ m/synonyms:/) &&
         ($Headers[$i] !~ m/gene_properties:/) &&
         ($Headers[$i] !~ m/snRNA symbol:/) &&
         ($Headers[$i] =~ m/pseudogene symbol:/))
       {
          $Headers[$i] =~ m/(C[A-Z][0-9]+) pseudogene symbol:([^ ]+) (FB[a-zA-Z0-9]+) /;
          $CG[$num] = $1;
          $ID[$num] = $2;
          $FB[$num] = $3;
#          print "$i\t\t$CG[$num] = $ID[$num] = $FB[$num] syn- prop- rna- pseudo+\n";
          $num++;

         UniqueList($1,$2,$3,$File_In,$Sequences[$i]);
       }
     elsif( (($Headers[$i] =~ m/FLAGGED_PROBLEMATIC/) ||
         ($Headers[$i] !~ m/synonyms:/) ||
         ($Headers[$i] !~ m/gene_properties:/) ) &&
         ($Headers[$i] !~ m/snRNA symbol:/) &&
         ($Headers[$i] =~ m/pseudogene symbol:/))
       {
          $Headers[$i] =~ m/(C[A-Z][0-9]+) pseudogene symbol:([^ ]+) .* (FB[a-zA-Z0-9]+) /;
          $CG[$num] = $1;
          $ID[$num] = $2;
          $FB[$num] = $3;
#          print "$i\t\t$CG[$num] = $ID[$num] = $FB[$num] syn- prop- snRNA- pseudo+ flagged+\n";
          $num++;

         UniqueList($1,$2,$3,$File_In,$Sequences[$i]);
       }
     elsif( ($Headers[$i] !~ m/FLAGGED_PROBLEMATIC/) &&
         ($Headers[$i] !~ m/synonyms:/) &&
         ($Headers[$i] !~ m/gene_properties:/) &&
         ($Headers[$i] =~ m/snRNA symbol:/))
       {
          $Headers[$i] =~ m/(CR[0-9]+) snRNA symbol:([^ ]+) (FB[a-zA-Z0-9]+) /;
          $CG[$num] = $1;
          $ID[$num] = $2;
          $FB[$num] = $3;
#          print "$i\t\t$CG[$num] = $ID[$num] = $FB[$num] syn- prop- rna+\n";
          $num++;

         UniqueList($1,$2,$3,$File_In,$Sequences[$i]);
       }
     elsif( (($Headers[$i] =~ m/FLAGGED_PROBLEMATIC/) ||
         ($Headers[$i] =~ m/synonyms:/) ||
         ($Headers[$i] =~ m/gene_properties:/) )&&
         ($Headers[$i] =~ m/snRNA symbol:/))
       {
          $Headers[$i] =~ m/(CR[0-9]+) snRNA symbol:([^ ]+) .* (FB[a-zA-Z0-9]+) /;
          $CG[$num] = $1;
          $ID[$num] = $2;
          $FB[$num] = $3;
#          print "$i\t\t$CG[$num] = $ID[$num] = $FB[$num] syn- prop- rna+\n";
          $num++;

         UniqueList($1,$2,$3,$File_In,$Sequences[$i]);
       }
     elsif( ($Headers[$i] !~ m/FLAGGED_PROBLEMATIC/) &&
         ($Headers[$i] !~ m/synonyms:/) &&
         ($Headers[$i] =~ m/gene symbol:/) &&
         ($Headers[$i] =~ m/gene_properties:/) )
       {
          $Headers[$i] =~ m/(CG[0-9]+) gene symbol:([^ ]+) .* (FB[a-zA-Z0-9]+) /;
          $CG[$num] = $1;
          $ID[$num] = $2;
          $FB[$num] = $3;
#          print "$i\t\t$CG[$num] = $ID[$num] = $FB[$num] syn- prop+\n";
          $num++;

         UniqueList($1,$2,$3,$File_In,$Sequences[$i]);
       }
     elsif( ($Headers[$i] !~ m/FLAGGED_PROBLEMATIC/) && 
            ($Headers[$i] =~ m/synonyms:/) &&
            ($Headers[$i] =~ m/gene symbol:/) )
       {
                                                                                          
          $Headers[$i] =~ m/(CG[0-9]+) gene symbol:([^ ]+) .* (FB[a-zA-Z0-9]+) /;
          $CG[$num] = $1;
          $ID[$num] = $2;
          $FB[$num] = $3;
#          print "$i\t\t$CG[$num] = $ID[$num] = $FB[$num] syn+\n";
          $num++;

         UniqueList($1,$2,$3,$File_In,$Sequences[$i]);
       }
     else
       {
          $Headers[$i] =~ m/(CG[0-9]+) gene symbol:([^ ]+).*FLAGGED_PROBLEMATIC (FB[a-zA-Z0-9]+) /;
          $CG[$num] = $1;
          $ID[$num] = $2;
          $FB[$num] = $3;
#          print "$i\t\t$CG[$num] = $ID[$num] = $FB[$num] PROBLEMATIC?\n";
          $num++;

         UniqueList($1,$2,$3,$File_In,$Sequences[$i]);
       }


     if( (length($CG[$num-1]) < 1) || (length($ID[$num-1]) < 1) || (length($FB[$num-1]) < 1))
       {
          print "error not picking up\n";
          printf "%s %d %s %d %s %d\n", $CG[$num], length($CG[$num]), 
                                        $ID[$num], length($ID[$num]), 
                                        $FB[$num], length($FB[$num]);
          print "$Headers[$i-1]\n\n$Headers[$i]\n\n";

          print "File_In = $File_In\n";

          exit;
       }
   }
undef(@CG);
undef(@ID);
undef(@FB);
   
}

#-------------------------------------------------------------------------
sub Deal_With_Input {
                                                                                                 
my($InputFile, @whole, $i, $header_count, $seq_count,  );
                                                                                                 
$InputFile = $_[0];

#print "+++$InputFile+++\n";
                                                                                                 
open(INPUT, "$InputFile") || die "Can't open $InputFile: $!\n";
@whole = <INPUT>;
close INPUT || die "Can't close $InputFile: $!\n";
                                                                                                 
for($i=0; $i<=$#whole; $i++)
   {
     #$whole[$i] = uc($whole[$i]);
                                                                                                 
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
     $Sequences[$i] =~ s/\n//g;
     $Sequences[$i] = uc($Sequences[$i]);
     $Headers[$i]   =~ s/\(GO.*//g;
     $Headers[$i]   =~ s/>//g;
   }
                                                                                                 
undef(@whole);
undef($InputFile);
undef($seq_count);
undef($header_count);
}
                                                                                                 
#-------------------------------------------------------------------------
# FIN
#-------------------------------------------------------------------------
