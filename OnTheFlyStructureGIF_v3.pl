#!/usr/bin/perl -w

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
use CGI qw(:standard);
use DBI;
use lib "/home/+user+/public_html/RNA/";
use Functions  qw(:ALL);

my( %HostStuff, $q, $IPs, $IP, $In_Phase, $In_UID, $In_Seq, $In_Size, 
    $In_Str, $In_Table, $In_Desc, );

$q = new CGI;

%HostStuff = HostSetUP($ENV{SERVER_ADDR});
$IPs       = $HostStuff{'IPs'};
$IP        = $HostStuff{'IP'};

$In_Phase = $q->param('Phase');  #ALGORITHM MFOLD OR LFOLD
$In_UID   = $q->param('UID');  #UID FOR LABEL
$In_Seq   = $q->param('Seq');    #SEQUENCE FOR DB LOOK UP
$In_Size  = $q->param('Size');   #IMAGE SIZE
$In_Str   = $q->param('Str');  
$In_Desc  = $q->param('Desc');   #TABLE NAME

#LogUser($0);

if( ($In_Phase =~ m/FOLD/) && ($In_UID !~ m/undefined/) )#&& ($In_Table =~ m/\b[MR]\b/) )
  {
     P_FromDB();
  }
elsif( ($In_Phase eq "RNAMOTIF") && ($In_UID !~ m/undefined/) )
  {
     P_FromDB();
#      PrintOutput_Error("$In_Str");
  }
else
   {
     PrintOutput_Error("No Structure Selected"); 
   }

#------------------------------------------------------------------------------
sub P_FromDB {

my( $i, @Seq, @Str, $dbh, $q, $SQL, $row, $q2, $SQL2, $row2, 
   $Style, $Size, $Output, $Sequence, $Structure, $Energy, $CG, $Header);

$Style = "bases";

$dbh = DBI->connect('DBI:mysql:RNASearch_r4_v1:+IP+:3306',
                    '+username+', '+password+',
                    { RaiseError => 1, AutoCommit => 1} )
                    || die "Database connection not made: $DBI::errstr";
               
$SQL2 = "SELECT CG, Sequence FROM $In_Desc WHERE UID='$In_UID' ";
$q2 = $dbh->prepare($SQL2);
$q2->execute;
while( $row2 = $q2->fetchrow_hashref() )
   {
     $CG = $row2->{'CG'};
     $CG =~ s/-.*//g;

     if( length($CG) > 12 )
       {
         $CG = substr($CG,0,9) . "...";
       }                                                                                       
     $Header   = "$CG UID=$In_UID $In_Phase $In_Desc";
     $Sequence = $row2->{'Sequence'};
   }

if($In_Phase =~ m/FOLD/)
  {
    $SQL = "SELECT $In_Phase\_SS, $In_Phase\_Energy FROM SecondaryStructures_Hits " . 
           "WHERE Sequence='$Sequence' "; 
    $q = $dbh->prepare($SQL);
    $q->execute;
    while( $row = $q->fetchrow_hashref() )
      {
        $Structure = $row->{"$In_Phase\_SS"};
        $Energy    = $row->{"$In_Phase\_Energy"};
        chomp($Energy);
        $Header   .= " E=$Energy";
      }
   }
elsif($In_Phase eq "RNAMOTIF")
  {
    $Sequence = "";
    @Seq    = split(//,$In_Seq);
    @Str    = split(//,$In_Str); 

    for($i=0; $i<=$#Str; $i++)
       {
         if($Str[$i] ne "-")
           {
             $Sequence  .= $Seq[$i];
             $Structure .= $Str[$i];
           }
       }

    $Energy = "-20.00";
    $Header = "$CG UID=$In_UID RNAMOTIF $In_Desc E=";
  }

$dbh->disconnect;

$Sequence = lc($Sequence);
open(NEWFASTA,">ss_image.fasta");
print NEWFASTA "$Sequence\n$Structure ($Energy)\n";
close NEWFASTA;

$ENV{'PATH'} = '/bin:/usr/bin:/usr/local/bin';

my(@Output);
if($In_Phase eq "RNAMOTIF")
  {
    @Output    = `RNAeval -T 25 < ss_image.fasta`;
    $Output[$#Output] =~ s/.* //g;
    $Output[$#Output] =~ s/[\(\)]//g;
    $Output[$#Output] =~ s/ //g;
    chomp($Output[$#Output]);
    $Header = $Header . $Output[$#Output];
  }


`cat ss_image.fasta | /usr/share/ViennaRNA/bin/b2ct > ss_image.ct`;
`cat /usr/local/lib/mfold-3.1/$Style.nav | sed 's/#.*//g' | sed 's/LAB_FR/20    /' | sed 's/ROT_ANG/0      /' | /usr/local/bin/naview.exe ss_image.ct ss_image.plt2`;
`sed 's/CTA     10.800 1.0.*/CTA     10.800 1.0 "$Header"/g' ss_image.plt2 > ss_image_2.plt2`; 
`/usr/local/bin/plt22gif -r $In_Size -o ss_image_M -x ss_image_2.plt2`;


$Output  = "<CENTER><TABLE STYLE='border:1px;border-style:solid' BGCOLOR=white WIDTH=55>";
$Output .= "<TR><TD>";
$Output .= "<FONT FACE=courier SIZE=1>>$Header<BR>";
$Output .= "<FONT COLOR=green>$Sequence</FONT><BR>";
$Output .= "<FONT COLOR=blue>$Structure</FONT><P>";

if($In_Phase =~ m/FOLD/)
  {
    $Output .= "Change Image Size: &nbsp";
    $Output .= "<A HREF='http://$IPs/cgi-rna/Version_3/OnTheFlyStructureGIF_v3.pl?" . 
               "Phase=MFOLD&Desc=$In_Desc&UID=$In_UID&Seq=$Sequence&Size=40' " . 
               "STYLE='TEXT-DECORATION: NONE'>[40]</A> ";
    $Output .= "<A HREF='http://$IPs/cgi-rna/Version_3/OnTheFlyStructureGIF_v3.pl?" .
               "Phase=MFOLD&Desc=$In_Desc&UID=$In_UID&Seq=$Sequence&Size=50' " . 
               "STYLE='TEXT-DECORATION: NONE'>[50]</A> ";
    $Output .= "<A HREF='http://$IPs/cgi-rna/Version_3/OnTheFlyStructureGIF_v3.pl?" .            
               "Phase=MFOLD&Desc=$In_Desc&UID=$In_UID&Seq=$Sequence&Size=72' " . 
               "STYLE='TEXT-DECORATION: NONE'>[72]</A> ";
    $Output .= "<A HREF='http://$IPs/cgi-rna/Version_3/OnTheFlyStructureGIF_v3.pl?" . 
               "Phase=MFOLD&Desc=$In_Desc&UID=$In_UID&Seq=$Sequence&Size=150' " .
               "STYLE='TEXT-DECORATION: NONE'>[150]</A> ";
  }
elsif($In_Phase eq "RNAMOTIF")
  {
    $Output .= "Change Image Size: &nbsp";
    $Output .= "<A HREF='http://$IPs/cgi-rna/Version_3/OnTheFlyStructureGIF_v3.pl?" .
               "Phase=RNAMOTIF&Desc=$In_Desc&UID=$In_UID&Seq=$Sequence&Str=$Structure&Size=40' " .
               "STYLE='TEXT-DECORATION: NONE'>[40]</A> ";
    $Output .= "<A HREF='http://$IPs/cgi-rna/Version_3/OnTheFlyStructureGIF_v3.pl?" .
               "Phase=RNAMOTIF&Desc=$In_Desc&UID=$In_UID&Seq=$Sequence&Str=$Structure&Size=50' " .
               "STYLE='TEXT-DECORATION: NONE'>[50]</A> ";
    $Output .= "<A HREF='http://$IPs/cgi-rna/Version_3/OnTheFlyStructureGIF_v3.pl?" .
               "Phase=RNAMOTIF&Desc=$In_Desc&UID=$In_UID&Seq=$Sequence&Str=$Structure&Size=72' " .
               "STYLE='TEXT-DECORATION: NONE'>[72]</A> ";
    $Output .= "<A HREF='http://$IPs/cgi-rna/Version_3/OnTheFlyStructureGIF_v3.pl?" .
               "Phase=RNAMOTIF&Desc=$In_Desc&UID=$In_UID&Seq=$Sequence&Str=$Structure&Size=150' " .
               "STYLE='TEXT-DECORATION: NONE'>[150]</A> ";
  }

$Output .= "</TD></TR><TR><TD><IMG SRC='http://$IP/Version_3/ss_image_M.gif'>";


$Output .= "</TD></TR></TABLE>";

print "Content-type: text/html\n\n";
                                                                                              
print <<END_OF_PAGE;
                                                                                              
<HTML>
<HEAD>
<TITLE>RNA Secondary Structure Image</TITLE>
</HEAD>
<BODY BGCOLOR=#9999FF onload="self.focus();">

<FONT FACE=courier COLOR=white SIZE=4>RNA Secondary Structure Image</FONT>
<BR> <BR>
                                                                                              
$Output
                                                                                              
</BODY>
</HTML>
END_OF_PAGE
}
