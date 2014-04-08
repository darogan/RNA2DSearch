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
use LWP::Simple;
use DBI;

my (@ParsedSeq, @BinaryImages, $i, @Matches, $CutOff, @RepSplit, $Counter,
    $dbh, $q1, $QuotedImage1, $QuotedImage2, $sql, @Annotation, 
    );

@Matches = ReadMatchesFile($ARGV[0]);

$ARGV[0] =~ m/(.*)_(.*)_dmel_RELEASE3-1.FASTA\..*ches_(.*)_(.*)_(.*)/;

print  "+----------+----------------------+-------------+--------+-------+\n";
print  "+ XSOMEARM + DATABASE             + TEMPERATURE + WINDOW + QUERY +\n";
print  "+----------+----------------------+-------------+--------+-------+\n";
printf "+ %8s + %20s + %11s + %6s + %6s + %7s +\n", $1, $2, $3, $4, $5;
print  "+----------+----------------------+-------------+--------+-------+\n";

my $XSome   = $1;
my $DB      = $2; # need to ensure that there are no underscores in database names!!!
my $Temp    = $3;
my $Window  = $4;
my $Query   = $5;

$Counter = 0;
$CutOff  = $ARGV[1];

$dbh = DBI->connect('DBI:mysql:RNALfold_Results:192.168.4.127:3306', '+username+', '+password+',
                    { RaiseError => 1, AutoCommit => 1} );
$dbh->{LongReadLen} = 16384;
$dbh->{LongTruncOk} = 1;


my($sql2, $q2, $random, $QuotedImage3,  );

my ($QuotedDB, $QuotedStart, $QuotedGene, $QuotedQuery, $QuotedLength, $QuotedScore,
    $QuotedEnergy, $QuotedXSome, $QuotedTemp, $QuotedWindow, $QuotedSequence, 
    @ScoreMultitude, @QuotedScoreMultitude, $QuotedAnnotation, $QuotedStructure, 
    $j, $Annotation );

for($i=0; $i<=$#Matches; $i++)
   {
     @RepSplit = split(/ +/, $Matches[$i][0]);

     $RepSplit[1] =~ s/ID=//g;  $RepSplit[1] =~ s/\(.*\)//g; 
     $RepSplit[2] =~ s/Length=//g;
     $RepSplit[3] =~ s/Start=//g;
     $RepSplit[4] =~ s/Score=//g; 

     #new to deal with the multitude of scoring parameters
     @ScoreMultitude = split(/\]\[/, $RepSplit[4] );

print "---$RepSplit[4]---\n";

     for($j=0; $j<=$#ScoreMultitude; $j++)
        {
          $ScoreMultitude[$j] =~ s/\[//g;
          $ScoreMultitude[$j] =~ s/\]//g;

          $QuotedScoreMultitude[$j] = $dbh->quote($ScoreMultitude[$j]);
        }

     $RepSplit[5] =~ s/Energy=//g; $RepSplit[5] =~ s/\[//g; $RepSplit[5] =~ s/\]//g;
    
     $Matches[$i][1] =~ s/SEQ: //g; 
     $Matches[$i][2] =~ s/SS_: //g;

     $QuotedGene       = $dbh->quote($RepSplit[1]);
     $QuotedStructure  = $dbh->quote($Matches[$i][2]);
     
     $QuotedXSome      = $dbh->quote($XSome);
     $QuotedDB         = $dbh->quote($DB);
     $QuotedQuery      = $dbh->quote($Query);
     $QuotedLength     = $dbh->quote($RepSplit[2]);
     $QuotedStart      = $dbh->quote($RepSplit[3]);
     $QuotedEnergy     = $dbh->quote($RepSplit[5]);
     $QuotedTemp       = $dbh->quote($Temp);
     $QuotedWindow     = $dbh->quote($Window);
     $QuotedSequence   = $dbh->quote($Matches[$i][1]);

     if($RepSplit[4] <= $CutOff)
       {
         $Annotation = "";
         @Annotation = GetGeneName($RepSplit[1]);
         for($j=0; $j<=$#Annotation; $j++)
            {
              $Annotation = $Annotation . $Annotation[$j] . "xxx";
            }
         $Annotation =~ s/xxx$//g;
         undef(@Annotation);

         print "$i - Annotation = $Annotation, $XSome\n";

         $QuotedAnnotation = $dbh->quote($Annotation);

         @BinaryImages = Mfolder($Matches[$i][1]);
         printf "There are %d images for sequence %d\n", 
                $#BinaryImages+1, $Counter+1;

         my $Bin = "None";

         if($#BinaryImages == 1)
           {
             $QuotedImage1 = $dbh->quote($BinaryImages[0]);
             $QuotedImage2 = $dbh->quote($BinaryImages[1]);
           }
         else
           {
             $QuotedImage1 = $dbh->quote($BinaryImages[0]);
             $QuotedImage2 = $dbh->quote($Bin);
           }

$sql = qq`INSERT INTO Matches_Score_All VALUES ('', $QuotedXSome, $QuotedDB, $QuotedGene, $QuotedAnnotation, $QuotedQuery, $QuotedLength, $QuotedStart,  $QuotedScore, $QuotedTemp, $QuotedWindow, $QuotedEnergy, $QuotedSequence, $QuotedStructure, $QuotedImage1, $QuotedImage2,0,'none') `;


#         $q1 = $dbh->prepare($sql);
#         $q1->execute;
#         $q1->finish;
         undef($q1);

         $Counter++;
       }

   }

$dbh->disconnect;

print "Within range = $Counter\n";

#----------------------------------------------------------------------------------------
sub ReadMatchesFile {

my ($MatchesFile, @MatchFile, $i, @Matches, $num, $num2, $num3, );

$MatchesFile = $_[0];

open (MATCHES, "$MatchesFile") || die "Can't open $MatchesFile :$!\n";
@MatchFile = <MATCHES>;
close MATCHES || die "Can't close $MatchesFile :$!\n";

for($i=0; $i<=$#MatchFile; $i++)
   {
      if($MatchFile[$i] =~ m/^(REP.*)/)
        {
          $Matches[$num][0] = $1;
          $num++;
        }
      elsif($MatchFile[$i] =~ m/^(SEQ.*)/)
        {
          $Matches[$num2][1] = $1;
          $num2++;
        }
      elsif($MatchFile[$i] =~ m/^(SS_.*)/)
        {
          $Matches[$num3][2] = $1;
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

#----------------------------------------------------------------------------------------
sub Mfolder {

my ($SeqFile, $rand, @Files, $i, $Output, $Count, $buffer, $binImage, 
    $ImgFile, $byteCount, @binImage, $ParsedSeq, );

srand;
$ENV{MFOLDLIB} = '/usr/local/lib/mfold-3.1';
$ENV{MFOLDBIN} = '/usr/local/bin';
$ParsedSeq     = $_[0];
$rand          = int(rand(5000));
$SeqFile       = "tmpfile_$rand";

#Convert sequence form DNA to RNA for mfold
$ParsedSeq =~ s/T/U/g;
$ParsedSeq =~ s/t/u/g;

open (SEQ, ">$SeqFile") || die "Can't open $SeqFile :$!\n";
print SEQ ">$SeqFile\n$ParsedSeq\n";
close SEQ || die "Can't close $SeqFile :$!\n";

$Output = `mfold SEQ=$SeqFile T=25 MAX=2 NA=RNA `;

@Files = <$SeqFile*>;

for($i=0; $i<=$#Files; $i++)
   {
     if($Files[$i] =~ m/$SeqFile.([0-9])\.gif\b/)
       {
          $Count   = $1;
          $ImgFile = $SeqFile . "_$Count.gif";

          open(IMG, "$ImgFile") || die "Can't open $ImgFile :$!\n";
          binmode(IMG);
          while(read(IMG,$buffer,1))
               {
                $binImage[$Count-1] = $binImage[$Count-1] . $buffer;
               }
          close IMG || die "Can't close $ImgFile :$!\n";
       }
     unlink($Files[$i]);
   }
return @binImage;
}

#------------------------------------------------------------------------------
sub GetGeneName {

my ($Contents, @Content, $i, $num, @Annotation, $AnnotFile, $GeneLink,  );
                                                                                
$GeneLink = $_[0];
                                                                                
$GeneLink =~ s/_.*//g;
$GeneLink =~ s/-.*//g;
                                                                                
print "+++$GeneLink+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

$AnnotFile = "/usr/people/+user+/Databases/DB_ANNOTATIONS/" . 
             $XSome .
             "_annotation_dmel_RELEASE3-1.FASTA.headers";

#open(ANNOT,"$AnnotFile") || die "Cant open $AnnotFile:$!\n";
#@GeneHeaders = <ANNOT>;
#close ANNOT || die "Cant close $AnnotFile:$!\n";

my $Annot = `grep ">$GeneLink" $AnnotFile`; 

if(length($Annot) < 1)
  {
    print "ANNOTATION MISSING FROM FILE\n";
    exit;
  }

if($GeneLink =~ m/CG/)
  {
    $Annot =~ m/gene symbol:([0-9A-Za-z\_:-]*) /;
    $Annotation[0] = $1;
    $Annot =~ m/(\(GO:.*\))/;
    $Annotation[1] = $1;
    $Annotation[1] =~ s/[\(\)]//g;
    $Annotation[1] =~ s/GO:[0-9]*//g;
    $Annotation[1] =~ s/" *"/,/g;
    $Annotation[1] =~ s/"$//g;
    $Annotation[1] =~ s/"//g;
    $Annotation[1] =~ s/^ *//g;
  }
elsif($GeneLink =~ m/TE/)
  {
    $Annot =~ m/transposable_element symbol:(.*) FBti/;
    $Annotation[0] = $1;
  }
elsif($GeneLink =~ m/CR/)
  {
    $Annot =~ m/symbol:([0-9A-Za-z\_:-]*) /;
    $Annotation[0] = $1;
    $Annot =~ m/(\(GO:.*\))/;
    $Annotation[1] = $1;
    $Annotation[1] =~ s/[\(\)]//g;
    $Annotation[1] =~ s/GO:[0-9]*//g;
    $Annotation[1] =~ s/" *"/,/g;
    $Annotation[1] =~ s/"$//g;
    $Annotation[1] =~ s/"//g;
    $Annotation[1] =~ s/^ *//g;
  }
else
  {
    print "ENTRY NOT CG TE or CG\n";
    exit;
  }

undef($Annot);

return @Annotation;
}

#------------------------------------------------------------------------------
sub GetGeneName_Old {

my ($Contents, @Content, $i, $num, @Annotation, $Link, $GeneLink,  );
                                                              
$GeneLink = $_[0];

$GeneLink =~ s/_.*//g;
$GeneLink =~ s/-.*//g;

print "+++$GeneLink\n";
                                                                                
if($GeneLink =~ m/CG/)
  {
    $GeneLink =~ s/CG//g;
    $Link = "http://flybase.bio.indiana.edu/.bin/fbidq.html?FBgn000$GeneLink"
  }
elsif($GeneLink =~ m/TE/)
  {
    $GeneLink =~ s/TE//g;
    $Link = "http://flybase.bio.indiana.edu/cgi-bin/fbidq.html?FBti00$GeneLink";
  }
elsif($GeneLink =~ m/CR/)
  {
    $GeneLink =~ s/CR//g;
    $Link = "http://flybase.bio.indiana.edu/.bin/fbidq.html?FBgn000$GeneLink";
  }
else
  {
     print "Unknown GeneLink ID ($GeneLink)\n";
     exit;
  }

                  
$Contents = LWP::Simple::get($Link);
@Content  = split(/\n/,$Contents);
                                                               
for($i=0; $i<=$#Content; $i++)
   {
     if($Content[$i] =~ m/Full name/)
       {
         $Content[$i] =~ s/<[^>]*>//g;
         $Content[$i] =~ s/.*Full name//g;
         $Content[$i] =~ s/^ *//;
         $Content[$i] =~ s/ *$//;
         $Annotation[$num] = $Content[$i];
         $num++;
       }
     elsif($Content[$i] =~ m/Annotation of/)
       {
         $Content[$i+7] =~ s/<[^>]*>//g;
         $Content[$i+7] =~ s/^ *//;
         $Content[$i+7] =~ s/ *$//;
         $Annotation[$num] = $Content[$i+7];
         $num++;
       }
     elsif($Content[$i] =~ m/TITLE/)
       {
         $Content[$i] =~ s/.*Report: //g;
         $Content[$i] =~ s/<[^>]*>//g;
         $Content[$i] =~ s/^ *//;
         $Content[$i] =~ s/ *$//;
         $Annotation[$num] = $Content[$i];
         $num++;
       }
   }
return @Annotation;
}

#----------------------------------------------------------------------------------------
# FIN
#----------------------------------------------------------------------------------------
