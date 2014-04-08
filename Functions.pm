package Functions;

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

# -------------------------------------------------------------------------
# Module for logging usage of RNA matches software
# Version 1
# -------------------------------------------------------------------------

use strict;
use warnings;
use DBI;

use vars qw( @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS 
             $VersionTable $ImageFileName $ImageMap %Clones
             $IPTest %HostDetails @RestrictionSites @Header @Tail 
             $MFOLD $LFOLD 
             %RNABinding %LocalizingGene @RNAMotifMatches);
use Exporter;
@ISA = qw(Exporter);

@EXPORT      = qw( VersionTable CreateMapImage QueryDB_Clones
                   QueryDB_RNABinding_Exists QueryDB_LocalisingGene_Exists
                   QueryDB_SecondaryStructures_Hits QueryDB_RNAMotif_r4_Descr_X 
                   HostSetUP MakeTail MakeHeader LogUser Authentification 
                   PrintOutput_Error FindRestrictionSites);
@EXPORT_OK   = qw( $VersionTable $ImageFileName $ImageMap %Clones
                   $IPTest %HostDetails @RestrictionSites @Header @Tail 
                   $MFOLD $LFOLD
                   @RNAMotifMatches %RNABinding %LocalizingGene);
%EXPORT_TAGS = ( ALL => [qw ( VersionTable 
                              CreateMapImage QueryDB_Clones
                              QueryDB_RNABinding_Exists
                              QueryDB_LocalisingGene_Exists
                              QueryDB_RNAMotif_r4_Descr_X
                              QueryDB_SecondaryStructures_Hits
                              HostSetUP LogUser Authentification 
                              MakeHeader MakeTail
                              PrintOutput_Error FindRestrictionSites)]
               );

#------------------------------------------------------------------------------
sub VersionTable {

my( $VersionTable, $i, @Tables, @TB_Entries, @TB_Desc, @TB_Notes, 
    $dbh, $q, @row, $row, $num, $SQL, $DB);

# DEALING WITH VERSIONS
if($_[0] eq 3){ $DB = "RNASearch_r4_v1";}
elsif($_[0] eq 2){ $DB = "RNASearch_r3_v4";}

$dbh = DBI->connect("DBI:mysql:$DB:+IP+:3306",
                    '+username+', '+password+',
                    { RaiseError => 1, AutoCommit => 1} )
                    || die "Database connection not made: $DBI::errstr";

# FIND OUT WHAT TABLES THERE ARE
$SQL = "SHOW TABLES LIKE 'RNAMotif\%Descr_%'";
$num = 0;
$q = $dbh->prepare($SQL);
$q->execute;
while( @row = $q->fetchrow)
 {
    $Tables[$num] = "@row";
    $num++;
 }
$q->finish;

# FIND OUT HOW POPULTAED THE TABLES ARE
for($i=0; $i<=$#Tables; $i++)
   {
     $SQL = "SELECT * FROM $Tables[$i]";
     $q = $dbh->prepare($SQL);
     $q->execute;

     $num = 0;
     while( $row = $q->fetchrow_hashref() )
          {
            $num++;
          }
     $TB_Entries[$i] = $num;
     $q->finish;


     $SQL = "SELECT Des, Notes from RNAMotif_TB_2_Desc WHERE TB = '$Tables[$i]' ";
     $q = $dbh->prepare($SQL);
     $q->execute;
     while( $row = $q->fetchrow_hashref() )
          {
            $TB_Desc[$i] = $row->{'Des'};
            $TB_Notes[$i] = $row->{'Notes'};
          }

     if( length($TB_Desc[$i]) < 1)
       {
         $TB_Desc[$i] = "<FONT COLOR=red>Description not updated yet!</FONT>";
       }

     $q->finish;
   }

$dbh->disconnect;

# PUT QUERIED DATA INTO A TABLE

$VersionTable  = "\n<TABLE STYLE='border-collapse:collapse;'>";
$VersionTable .= "<TR><TD ALIGN=center><FONT SIZE=1 FACE='sans,arial'><B>Description</B></FONT></TD>" .
                 "<TD ALIGN=center><FONT SIZE=1 FACE='sans,arial'><B>Entries</B></FONT></TD>" .
                 "<TD ALIGN=center><FONT SIZE=1 FACE='sans,arial'><B>DB ID</B></FONT></TD>" . 
                 "<TD ALIGN=center WIDTH=250><FONT SIZE=1 FACE='sans,arial'><B>Notes</B></FONT></TD>" .
                 "</TR>";


for($i=0; $i<=$#Tables; $i++)
   {
     $VersionTable .= "<TR>" . 
                      "<TD style='border:1px solid;padding-left:5px;padding-right:5px' VALIGN=top>" . 
                        "<FONT SIZE=1 FACE='sans,arial'>$TB_Desc[$i]</FONT></TD>" . 
                      "<TD style='border:1px solid;padding-left:5px;padding-right:5px' " . 
                        "VALIGN=top ALIGN=right>" . 
                        "<FONT SIZE=1 FACE='sans,arial' COLOR=blue>$TB_Entries[$i]</FONT></TD>" . 
                      "<TD style='border:1px solid;padding-left:5px;padding-right:5px' VALIGN=top>" . 
                        "<FONT SIZE=1 FACE='sans,arial'>$Tables[$i]</FONT></TD>" . 
                      "<TD style='border:1px solid;padding-left:5px;padding-right:5px' " . 
                        "VALIGN=top WIDTH=250>" . 
                        "<FONT SIZE=1 FACE='sans,arial'>$TB_Notes[$i]</FONT></TD>" .
                      "</TR>";
   }
$VersionTable .= "</TABLE>\n";

return $VersionTable;
}

#------------------------------------------------------------------------------
sub CreateMapImage {

my( $ImageFileName, $Table, @Table, $i, $j, @Bits, %Structure, $AAsPixel,
    @Decision, %Decision, $Height, $Width, $key, %Colours, $im, $LineCount,
    @Split, @OrderToPlot, %ColourChart, $Title, $ImageMap, $Coords,
  );

use GD;

$Table    = $_[0];
@Table    = split(/\n/,$Table);
$Title    = $_[1];
@Decision = split(/\n/,$_[2]);

for($i=0; $i<=$#Table; $i++)
   {
     $Table[$i]     =~ m/(.*)\t(.*)/;
     $Structure{$1} = $2;
   }

for($i=0; $i<=$#Decision; $i++)
   {
     $Decision[$i] =~ m/(.*)\t(.*)/;
     $Decision{$1} = $2;
   }

$Structure{'GN'} =~ s/[\[\]\-0]//g;
$AAsPixel        =  500 / $Structure{'GN'};

$Height = ($#Table * 20 ) + 100;
$Width  = 675;

$ImageMap = "<MAP NAME='GeneMapping'>";
$im = new GD::Image($Width,$Height);

$Colours{'white'}  = $im->colorAllocate(255,255,255);
$Colours{'black'}  = $im->colorAllocate(0,0,0);
$Colours{'red'}    = $im->colorAllocate(255,0,0);
$Colours{'blue'}   = $im->colorAllocate(0,0,255);
$Colours{'cyan'}   = $im->colorAllocate(0,255,255);
$Colours{'green'}  = $im->colorAllocate(0,255,0);
$Colours{'purple'} = $im->colorAllocate(200,0,200);
$Colours{'orange'} = $im->colorAllocate(255,153,0);

$im->transparent($Colours{'white'});
$im->interlaced('true');

$im->string(gdSmallFont,5,5,$Title,$Colours{'black'});

$LineCount = 1;
$im->string(gdSmallFont,15,($LineCount*25)+2.5,"Genomic",$Colours{'blue'});
$im->filledRectangle(125,$LineCount*25,125+($AAsPixel*$Structure{'GN'}),($LineCount*25)+10,$Colours{'blue'});
$LineCount++;

@OrderToPlot = ("CLONE", "TS", "CDS", "MATCH");
%ColourChart = ("CLONE", "cyan", "TS", "green",  "CDS", "green", "MATCH", "red");



for($j=0; $j<=$#OrderToPlot; $j++)
   {
     for $key (sort keys %Structure)
        {
          if($key =~ m/^$OrderToPlot[$j]/)
            {
              $im->filledRectangle(125,($LineCount*25)+4.5,
                                   625,($LineCount*25)+5.5,
                                        $Colours{'black'});
              $im->string(gdSmallFont,15,($LineCount*25),"$key",
                          $Colours{$ColourChart{$OrderToPlot[$j]}});
              @Split = split(/\s+|\t/,$Structure{$key});
              for($i=0; $i<=$#Split; $i++)
                 {
                   $Split[$i] =~ m/\[(.*)-(.*)\]/;
                   $im->filledRectangle(125+($AAsPixel*$1),$LineCount*25,
                                        125+($AAsPixel*$2),($LineCount*25)+10,
                                        $Colours{$ColourChart{$OrderToPlot[$j]}});
                   $Coords = sprintf("%.0f,%.0f,%.0f,%.0f",
                             125+($AAsPixel*$1),$LineCount*25,
                             125+($AAsPixel*$2),($LineCount*25)+10);

                   $Title = "$key [Range=$1-$2] [Exons=" . ($#Split+1) . "]";

                   $ImageMap .= "\n<AREA SHAPE='rect' HREF='' TITLE='$Title' " .
                                "onMouseover='Myform.Report.value=\"$Title $Decision{$key}\";" .
                                "return false;' " .
                                "COORDS='$Coords'>\n";
                 }
              $LineCount++;
              undef(@Split);
            }
        }
  }

$im->string(gdTinyFont,5,($LineCount*25),"Scale: 100nt", $Colours{'black'});
$im->filledRectangle(40,($LineCount*25)+10,
                     40+($AAsPixel*100),($LineCount*25)+11,
                     $Colours{'black'});

$im->string(gdTinyFont,500,($LineCount*25)+2.5,"Copyright R. S. Hamilton 2005", $Colours{'black'});

$ImageMap .= "</MAP><P>";

$ImageFileName = "NewAlignmentImage.png";
open(IMG, ">$ImageFileName") or die $!;
  binmode IMG;
print IMG $im->png;
close(IMG);

return ($ImageFileName, $ImageMap);
}

#-----------------------------------------------------------------------
sub QueryDB_SecondaryStructures_Hits {

my( $Sequence, $dbh, $SQL, $q, $row, $MFOLD, $LFOLD, $MEnergy, $LEnergy );

$Sequence = $_[0];

$dbh = DBI->connect('DBI:mysql:RNASearch_r4_v1:+IP+:3306',
                    '+username+', '+password+',
                    { RaiseError => 1, AutoCommit => 1} )
                    || die "Database connection not made: $DBI::errstr";

$SQL = "SELECT MFOLD_SS, LFOLD_SS, MFOLD_Energy, LFOLD_Energy " . 
       "FROM SecondaryStructures_Hits " . 
       "WHERE Sequence='$Sequence'";

$q = $dbh->prepare("$SQL");
$q->execute;
while( $row = $q->fetchrow_hashref() )
     {
       $MFOLD   = $row->{'MFOLD_SS'};
       $LFOLD   = $row->{'LFOLD_SS'};
       $MEnergy = $row->{'MFOLD_Energy'};
       $LEnergy = $row->{'LFOLD_Energy'};
     }
$q->finish;

$dbh->disconnect;
return ($MFOLD,$LFOLD,$MEnergy,$LEnergy);
}

#-----------------------------------------------------------------------
sub QueryDB_RNAMotif_r4_Descr_X {

my(@RNAMotifMatches, $SQL, $dbh, $q, $row, $num,  );

$SQL = $_[0];

$dbh = DBI->connect('DBI:mysql:RNASearch_r4_v1:+IP+:3306',
                    '+username+', '+password+',
                    { RaiseError => 1, AutoCommit => 1} )
                    || die "Database connection not made: $DBI::errstr";

$q = $dbh->prepare("$SQL");
$q->execute;
while( $row = $q->fetchrow_hashref() )
 {
   $RNAMotifMatches[$num]{'UID'}             = $row->{'UID'};
   $RNAMotifMatches[$num]{'Sequence'}        = $row->{'Sequence'};
   $RNAMotifMatches[$num]{'CG'}              = $row->{'CG'};
   $RNAMotifMatches[$num]{'Chromosome'}      = $row->{'Chromosome'};
   $RNAMotifMatches[$num]{'Annotation'}      = $row->{'Annotation'};
   $RNAMotifMatches[$num]{'Scores'}          = $row->{'Scores'};
   $RNAMotifMatches[$num]{'RNAMotif_Score'}  = $row->{'RNAMotif_Score'};
   $RNAMotifMatches[$num]{'RNAMotif'}        = $row->{'RNAMotif'};    
   $RNAMotifMatches[$num]{'RNAMotif_SS'}     = $row->{'RNAMotif_SS'};
   $RNAMotifMatches[$num]{'Location'}        = $row->{'Location'};

   $RNAMotifMatches[$num]{'Complexity'}      = $row->{'Complexity'};
   $RNAMotifMatches[$num]{'OldCandidate'}    = $row->{'OldCandidate'};

   $RNAMotifMatches[$num]{'CDS'}             = $row->{'CDS'};
   $RNAMotifMatches[$num]{'five_prime_UTR'}  = $row->{'five_prime_UTR'};    
   $RNAMotifMatches[$num]{'gene'}            = $row->{'gene'};
   $RNAMotifMatches[$num]{'intron'}          = $row->{'intron'};
   $RNAMotifMatches[$num]{'miscRNA'}         = $row->{'miscRNA'};
   $RNAMotifMatches[$num]{'pseudogene'}      = $row->{'pseudogene'};
   $RNAMotifMatches[$num]{'three_prome_UTR'} = $row->{'three_prome_UTR'};    
   $RNAMotifMatches[$num]{'transcript'}      = $row->{'transcript'};
   $RNAMotifMatches[$num]{'transposon'}      = $row->{'transposon'};

   $RNAMotifMatches[$num]{'Candidate'}       = $row->{'Candidate'};
   $RNAMotifMatches[$num]{'Comments'}        = $row->{'Comments'};    
   $RNAMotifMatches[$num]{'ByEye'}           = $row->{'ByEye'};

   if($RNAMotifMatches[$num]{'OldCandidate'} == 0)
     { $RNAMotifMatches[$num]{'OldCandidate'} = "No"}
   elsif($RNAMotifMatches[$num]{'OldCandidate'} == 1)
     { $RNAMotifMatches[$num]{'OldCandidate'} = "Gene Only"}
   elsif($RNAMotifMatches[$num]{'OldCandidate'} == 2)
     { $RNAMotifMatches[$num]{'OldCandidate'} = "Seq Only"}
   elsif($RNAMotifMatches[$num]{'OldCandidate'} == 3)
     { $RNAMotifMatches[$num]{'OldCandidate'} = "Gene + Seq"}
    $num++;
 }
$q->finish;
$dbh->disconnect;

return @RNAMotifMatches;
}

#-----------------------------------------------------------------------
sub QueryDB_LocalisingGene_Exists {

my(%LocalizingGene, $dbh, $q, $row,);

$dbh = DBI->connect('DBI:mysql:RNASearch_General:+IP+:3306',
                    '+username+', '+password+',
                    { RaiseError => 1, AutoCommit => 1} )
                    || die "Database connection not made: $DBI::errstr";

$q = $dbh->prepare("SELECT CG FROM LocalizationGene");
$q->execute;
while( $row = $q->fetchrow_hashref() )
 {
   $LocalizingGene{"$row->{'CG'}"} = 1;
 }
$q->finish;
$dbh->disconnect;

undef($dbh);
undef($q);
undef($row);

return %LocalizingGene;
}

#-----------------------------------------------------------------------
sub QueryDB_Clones {

my(%Clones, $dbh, $q, $row, );

$dbh = DBI->connect('DBI:mysql:RNASearch_General:+IP+:3306',
                    '+username+', '+password+',
                    { RaiseError => 1, AutoCommit => 1} )
                    || die "Database connection not made: $DBI::errstr";

$q = $dbh->prepare("SELECT CID, Sequence FROM Clones_Gold");
$q->execute;
while( $row = $q->fetchrow_hashref() )
 {
   $Clones{$row->{'CID'}} = $row->{'Sequence'};
 }
$q->finish;
$dbh->disconnect;

undef($dbh);
undef($q);
undef($row);

return %Clones;
}

#-----------------------------------------------------------------------
sub QueryDB_RNABinding_Exists {

my(%RNABinding, $dbh, $q, $row, );

$dbh = DBI->connect('DBI:mysql:RNASearch_General:+IP+:3306',
                    '+username+', '+password+',
                    { RaiseError => 1, AutoCommit => 1} )
                    || die "Database connection not made: $DBI::errstr";

$q = $dbh->prepare("SELECT CG FROM RNABinding");
$q->execute;
while( $row = $q->fetchrow_hashref() )
 {
   $RNABinding{"$row->{'CG'}"} = 1;
 }
$q->finish;
$dbh->disconnect;

undef($dbh);
undef($q);
undef($row);

return %RNABinding;
}

#-----------------------------------------------------------------------
sub MakeHeader {

my($IPs, $IP, @Header, $i, );

$IPs = $_[0];
$IP  = $_[1];

if($IPs eq "+IP+")
  {
    open(HEADER,"header.text");
  }
elsif($IPs eq "129.xxx.xxx.xxx")
  {
    open(HEADER,"header.text");
  }

@Header = <HEADER>;
close HEADER; 

for($i=0; $i<=$#Header; $i++)
   {      
     $Header[$i] =~ s/\bIPs\b/$IPs/g;
     $Header[$i] =~ s/\bIP\b/$IP/g;
   }
return @Header;
}

#-----------------------------------------------------------------------
sub MakeTail {

my(@Tail, $i, $IPs, $IP,);

$IPs = $_[0];
$IP  = $_[1];

if($IPs eq "+IP+")
  {
    open(TAIL,"tail.text");
  }
elsif($IPs eq "+IP+")
  {
    open(TAIL,"/home/ilan/public_html/RNASearch/General/tail.text");
  }

@Tail = <TAIL>;
close TAIL;

for($i=0; $i<=$#Tail; $i++)
   {
     $Tail[$i] =~ s/\bIPs\b/$IPs/g;
     $Tail[$i] =~ s/\bIP\b/$IP/g;
   }

return @Tail;
}

#-----------------------------------------------------------------------
sub FindRestrictionSites {

my ($Sequence, @RestrictionSites, $dbh, $q, $row, @Results, $num, 
    $matchStart, $instrumentedPattern, $nextStart, $Pattern, 
    @RestrictionSites_S, $i, $count,  );

$Sequence = $_[0];
use re qw(eval);
$dbh = DBI->connect('DBI:mysql:RNASearch_General:+IP+:3306',
                    '+username+', '+password+',
                    { RaiseError => 1, AutoCommit => 1} )
                    || die "Database connection not made: $DBI::errstr";

$q = $dbh->prepare("SELECT * FROM RestrictionEnzymes WHERE Active = '1'");
$q->execute;
$num=0;
while( $row = $q->fetchrow_hashref() )
     {
        $Results[$num][0] = $row->{'EnzymeName'};
        $Results[$num][1] = $row->{'Prototype'};
        $Results[$num][2] = $row->{'RecognitionSite'};
        $Results[$num][3] = $row->{'MethylationSite'};
        $Results[$num][4] = $row->{'PatternMatch'};
        $Results[$num][5] = $row->{'Active'};
        $num++;
     }
$q->finish;
$dbh->disconnect;

for($i=0; $i<=$#Results; $i++)
   {
     local $_;
     $_ = $Sequence;

     $Pattern = $Results[$i][4];
     $instrumentedPattern = qr/(?{ $matchStart = pos() })$Pattern/;

     while(/$instrumentedPattern/g)
          {
            $nextStart = pos()-1;

            push @RestrictionSites, { Enzyme       => $Results[$i][0], 
                                      CleavageSite => $Results[$i][2], 
                                      Start        => $matchStart,
                                      Stop         => $nextStart        };
          }                                                                                           
     undef($_);
     undef($nextStart);
     undef($matchStart);
     undef($instrumentedPattern);
   }

@RestrictionSites_S = sort { $a-> {Start} <=> $b->{Start} } @RestrictionSites;

undef($Sequence);
undef(@Results);
undef($num);

return @RestrictionSites_S;
}

#-----------------------------------------------------------------------
sub HostSetUP {

my(%HostDetails, $SERVER_ADDR, );

$SERVER_ADDR = $_[0];

if( $SERVER_ADDR eq "+IP+" )
  {
    $HostDetails{'IPs'}       = "+IP+";
    $HostDetails{'IP'}        = "+IP+/~+user+/RNA";
    $HostDetails{'BLAST_Loc'} = "/usr/progs/blast-2.2.11/bin";
    $HostDetails{'GD_Ext'}    = "png";
  }
elsif( $SERVER_ADDR eq "127.0.0.1" )
  {
    $HostDetails{'IPs'}       = "127.0.0.1";
    $HostDetails{'IP'}        = "127.0.0.1/~+user+/RNA";
    $HostDetails{'BLAST_Loc'} = "/usr/progs/BLAST";
    $HostDetails{'GD_Ext'}    = "png";
  }

return %HostDetails;
}

#-----------------------------------------------------------------------
sub LogUser {

my($dbh, $q, $q2, $q3, $q4, $row, $date, $num, $UserIP, $Program, $Program_D, 
   $Program_C, $UsersCount,   );

$UserIP     =  $ENV{'REMOTE_ADDR'};
$UsersCount =  0;
$Program    =  $_[0];
$Program    =~ s/.*\///g;
$Program    =~ s/.pl//;
$Program_D  =  $Program . "_D";
$Program_C  =  $Program . "_C";

$dbh = DBI->connect('DBI:mysql:RNALfold_Results:+IP+:3306',
                    '+username+', '+password+',
                    { RaiseError => 1, AutoCommit => 1} )
                    || die "Database connection not made: $DBI::errstr";

$q = $dbh->prepare("SELECT $Program_C, $Program_D FROM UserCounters WHERE UIP = '$UserIP'");
$q->execute;
while( $row = $q->fetchrow_hashref() )
  {
    $UsersCount = $row->{"$Program_C"};
  }
$q->finish;

if( $UsersCount eq undef)
  {
    $q2 = $dbh->prepare("INSERT INTO UserCounters VALUES ('$UserIP','','','','','','','','','','','','','','','','')");
    $q2->execute;
    $q2->finish;
    $UsersCount = 0;
  }

$UsersCount = $UsersCount + 1;
$q3 = $dbh->prepare("UPDATE UserCounters set $Program_C='$UsersCount' where UIP='$UserIP'");
$q3->execute;
$q3->finish;

$q4 = $dbh->prepare("UPDATE UserCounters set $Program_D=now() where UIP='$UserIP'");
$q4->execute;
$q4->finish;

$dbh->disconnect;
}

#------------------------------------------------------------------------------
sub Authentification {
my( $Output1, $Output2, @keys, $IPAddress, $num, $key, @user_info,
    @user_info_final, @IPSplit, $IPTest, $i,  );
                                                                                       
@keys = sort(keys(%ENV));
$num=0;
foreach $key (@keys)
   {
     $user_info[$num] = sprintf "$key = $ENV{$key}<BR>";
                                                                                        
                                                                                        
     if($key =~ /\bREMOTE_ADDR\b/)
       {
         $user_info_final[0] = $ENV{$key};
       }
     $num++;
   }
                                                                                        
@IPSplit = split(/\./, $user_info_final[0]);
                                                                                        
if(($IPSplit[0] == 127) && ($IPSplit[1] == 0))
  {
   $IPTest = "Pass";
  }
else
  {
    $IPTest = "Pass";
  }
                                                                                        
$IPAddress = $user_info_final[0];
return $IPTest;
}

#------------------------------------------------------------------------------
sub PrintOutput_Error {

my $Message = $_[0];
print "Content-type: text/html\n\n";

print <<END_OF_PAGE;                                                                                                    
<HTML>
<HEAD>
<TITLE>RNAMotif Search Results: Russell S. Hamilton</TITLE>
</HEAD>
<BODY BGCOLOR=#9999FF><FONT FACE="helvetica" SIZE=2>
<FONT FACE=courier COLOR=white SIZE=5>RNAMotif Search Results</FONT><P>
<TABLE STYLE='border:1px;border-style:solid' BGCOLOR=white>
<TR>
   <TD WIDTH=50 VALIGN=top><FONT FACE="helvetica" SIZE=3>Error:</TD>
   <TD VALIGN=top><FONT FACE="helvetica" SIZE=3>$Message<BR>
   <FONT COLOR=red>Access Denied</FONT></TD>
</TR>
</TABLE>
</BODY>
</HTML>
END_OF_PAGE

exit;
}

#--------------------------------------------------------------------------
1; # returns true for the use calls
# End of Program
#--------------------------------------------------------------------------
