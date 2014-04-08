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
use CGI qw(:standard);
use DBI;
use lib "/home/+user+/public_html/RNA";
use Functions  qw(:ALL);
                                                                                        
                                                                                        
my($i, @Header, @Tail, $IPs, $IP, $q, $IP_Test, $Results, $SQL, $Table, 
   %HostStuff, $In_LimitLow, $In_LimitHigh, $In_Form, $Range, 
   $In_QType, $In_RNAMotifDesc, $In_Filter, $In_FilterCO, $In_RNAMotif_Q, 
   $In_SeqIdentity, $In_Location, 
 );
    
%HostStuff = HostSetUP($ENV{SERVER_ADDR});
$IPs       = $HostStuff{'IPs'};
$IP        = $HostStuff{'IP'};
                                                                                    
$q = new CGI;
                                                                                        
@Header = MakeHeader($IPs,$IP);
@Tail   = MakeTail($IPs,$IP);                                                                
 
$IP_Test       = Authentification();
if($IP_Test =~ m/Pass/)
  {
    LogUser("RNAMotifResults");
    $In_Form         = $q->param('FormID');
    $In_QType        = $q->param('QType');
    $In_RNAMotifDesc = $q->param('RNAMotifDesc');
    $In_Filter       = $q->param('Filter');
    $In_FilterCO     = $q->param('FilterCutOff');
    $In_LimitLow     = $q->param('LimitLow');
    $In_LimitHigh    = $q->param('LimitHigh');
    $In_RNAMotif_Q   = $q->param('RNAMotif_Q');
    $In_SeqIdentity  = $q->param('SeqIdentity');
    $In_Location     = $q->param('Location');

    ($SQL,$Range) = DetermineSearchType($In_QType,$In_RNAMotifDesc,$In_LimitLow,$In_LimitHigh,
                                        $In_RNAMotif_Q,$In_Filter, $In_FilterCO);
    $Results = QueryTheDatabase($SQL,$Range);
    PrintOutput_Pass($Results);    
  }
else
  {
     PrintOutput_Error("Access Denied");
  }


#------------------------------------------------------------------------------
sub DetermineSearchType {

my($SQL, $Range, $QType, $RNAMotifDesc, $LimitLow, $LimitHigh,  
   $Filter, $FilterCO, $RNAMotif_Q, $Complexity, );

$QType        = $_[0];
$RNAMotifDesc = $_[1];
$LimitLow     = $_[2];
$LimitHigh    = $_[3];
$RNAMotif_Q   = $_[4];
$Filter       = $_[5];
$FilterCO     = $_[6];


if( $LimitLow >= $LimitHigh)
  {
    PrintOutput_Error("Invalid Range Entered");
  }
if( ($LimitHigh - $LimitLow) > 199)
  {
    PrintOutput_Error("Invalid Range Entered - too many!");
  }

if($Filter eq "LC")
  {
    $Complexity = " AND Complexity <= '$FilterCO' ";
  }
else
  {
    $Complexity = "";
  }

if($RNAMotifDesc eq "GandIwScore.rnamotif.desc")
  {
    $Table = "RNAMotif_r4_Descr_A";
  }
elsif($RNAMotifDesc eq "GandIwScoreUU.rnamotif.desc")
  {
    $Table = "RNAMotif_r4_Descr_B";
  }
elsif($RNAMotifDesc eq "GandIwScoreUU_CAA.rnamotif.desc")
  {
    $Table = "RNAMotif_r4_Descr_C";
  }
elsif($RNAMotifDesc eq "GonlyFixed.rnamotif.desc")
  {
    $Table = "RNAMotif_r4_Descr_D";
  }
elsif($RNAMotifDesc eq "IonlyFixed.rnamotif.desc")
  {
    $Table = "RNAMotif_r4_Descr_E";
  }
elsif($RNAMotifDesc eq "JockeyandIwScore.rnamotif.desc")
  {
    $Table = "RNAMotif_r4_Descr_F";
  }
elsif($RNAMotifDesc eq "K10andOrbwScore.desc")
  {
    $Table = "RNAMotif_r4_Descr_G";
  }
elsif($RNAMotifDesc eq "K10andOrbwScoreSeqCons.desc")
  {
    $Table = "RNAMotif_r4_Descr_H";
  }


my $Limiter;
if($LimitLow == 0)
  {
    $Limiter = "$LimitLow,$LimitHigh";
  }
else
  {
    $Limiter = sprintf("%d,%d", $LimitLow-1, $LimitHigh-$LimitLow);
  } 

if($QType eq "All")
  {
     $SQL = "SELECT * FROM $Table WHERE UID > 0 $Complexity ORDER BY CG, Sequence LIMIT $Limiter"; 
  }
elsif($QType eq "CG")
  {
     $SQL = "SELECT * FROM $Table WHERE CG LIKE 'CG%' $Complexity ORDER BY CG, Sequence LIMIT $Limiter";
  }
elsif($QType eq "TE")
  {
    $SQL = "SELECT * FROM $Table WHERE CG LIKE 'TE%' $Complexity ORDER BY CG, Sequence LIMIT $Limiter";
  }
elsif($QType eq "CR")
  {
    $SQL = "SELECT * FROM $Table WHERE CG LIKE 'CR%' $Complexity ORDER BY CG, Sequence LIMIT $Limiter";
  }
elsif($QType eq "IndividualID")
  {
    if( ($RNAMotif_Q !~ m/^[CT][GETR]/) || ($RNAMotif_Q =~ m/[^0-9a-zA-Z]/) )
      {
        PrintOutput_Error("Invalid Accession Number"); 
      }
    else
      {
        $SQL = "SELECT * FROM $Table WHERE CG LIKE '$RNAMotif_Q' " . 
               "$Complexity ORDER BY CG, Sequence LIMIT $Limiter";
      }
  }
elsif($QType eq "IndividualName")
  {
    if($RNAMotif_Q =~ m/[^0-9a-zA-Z{}]/)
      {
        PrintOutput_Error("Invalid Name: Contains illegal Characters");
      }
    else
      {
	$RNAMotif_Q = "%" . $RNAMotif_Q . "%";
        $SQL = "SELECT * FROM $Table WHERE Annotation LIKE '$RNAMotif_Q' " .                
               "$Complexity ORDER BY CG, Sequence LIMIT $Limiter";
      }
  }


$Range = "$LimitLow to $LimitHigh";
return $SQL, $Range;
}

#------------------------------------------------------------------------------
sub QueryTheDatabase {

my($Range, $Output, $i, $Rank, $SQL, @RNAMotifMatches, $dbh, $num, 
   $SeqPanel, %RNABinding, %Localizing, @DistanceScores,  
   $RNAMotifCount, $j, $RNAMotifDropDown, @RNAMotif_SSs, $MEnergy, $LEnergy, 
   @Comments, @ByEye, $CandidateString, $Candidate, $row, $q, 
   $O_Mid, $LinksImages, $UpdateBox, $O_Left, $O_Top, $MFOLD, $LFOLD, );

$SQL      = $_[0];
$Range    = $_[1];

$dbh = DBI->connect('DBI:mysql:RNASearch_r4_v1:+IP+:3306',
                    '+username+', '+password+',
                    { RaiseError => 1, AutoCommit => 1} )
                    || die "Database connection not made: $DBI::errstr";

@RNAMotifMatches = QueryDB_RNAMotif_r4_Descr_X($SQL);
%RNABinding      = QueryDB_RNABinding_Exists();
%Localizing      = QueryDB_LocalisingGene_Exists();

$Output = sprintf("<BR><FONT FACE=courier SIZE=3 COLOR=white>Range: %s Total displayed = %d<BR>",
          $Range, $#RNAMotifMatches+1);

$Output .= "<FONT FACE=courier SIZE=2>[grk +ve score][I factor -ve score]</FONT>\n<P>\n";

$Output .= "+++$SQL+++$#RNAMotifMatches+++\n<P>\n";

for($i=0; $i<=$#RNAMotifMatches; $i++)
   {
     ($MFOLD, $LFOLD, $MEnergy, $LEnergy) = QueryDB_SecondaryStructures_Hits($RNAMotifMatches[$i]{'Sequence'});

     if(length($MFOLD) < 1) {$MFOLD = "Coming Soon...";}
     if(length($LFOLD) < 1) {$LFOLD = "Coming Soon...";}

     $Rank = $i+1;
     $RNAMotifMatches[$i]{'RNAMotif'}   =~ s/xxx/\n/g;
     $RNAMotifMatches[$i]{'RNAMotif'}   =~ s/tmp .*[0-9] //g;
     $RNAMotifMatches[$i]{'RNAMotif'}   =  uc($RNAMotifMatches[$i]{'RNAMotif'});
     $RNAMotifMatches[$i]{'Annotation'} =~ s/xxx/ : /g;

     $RNAMotifMatches[$i]{'Scores'} =~ s/^\[//;
     @DistanceScores = split(/\]\[/,$RNAMotifMatches[$i]{'Scores'});


     ### START UPDATE BOX ###
     $Candidate    = $RNAMotifMatches[$i]{'Candidate'};
     $ByEye[$i]    = $RNAMotifMatches[$i]{'ByEye'};
     $Comments[$i] = $RNAMotifMatches[$i]{'Comments'};

     if($Candidate == 1) { $CandidateString = "CHECKED=checked"; }
     else { $CandidateString = ""; }

     $UpdateBox = "<TABLE CELLPADDING=1 CELLSPACING=1 " .
                  "STYLE='border:1px;border-style:solid;background-color:white'>" .

                  "\n<TR><TD ROWSPAN=2 ALIGN=center>" .
                  "<TEXTAREA NAME='com$i' ROWS=2 COLS=18 " .
                  "STYLE='font-family:sans,arial;font-size:6pt;color:blue;background-color:gold' " .
                  #" OnChange=\"UpdateComments(com$i,$RNAMotifMatches[$i]{'UID'});\" " .
                  ">$Comments[$i]</TEXTAREA></TD>" .

                  "\n<TD VALIGN=middle ALIGN=left>" .
                  "<IMG SRC='http://$IP/Images/syringe_resize.gif'></TD>" .

                  "<TD ALIGN=left><INPUT TYPE='CHECKBOX' NAME='can$i' $CandidateString " .
                  "STYLE='font-family:courier;font-size:7pt;color:blue;background-color:gold' " .
                  #"OnChange=\"UpdateCandidate(can$i,$RNAMotifMatches[$i]{'UID'});\" " .
                  "></TD>" . 

                  "<TD STYLE='border:1px;border-style:solid;background-color:white'>" . 
                  "<FONT FACE=sans,arial SIZE=1>Old Candidate: " . 
                  "<FONT COLOR=blue>$RNAMotifMatches[$i]{'OldCandidate'}</FONT></FONT></TD></TR>" .

                  "\n<TR><TD ALIGN=left VALIGN=middle><FONT FACE=courier SIZE=2>" .
                  "<IMG SRC='http://$IP/Images/eye.gif'></TD>" .

                  "<TD ALIGN=center><FONT FACE=courier SIZE=2 COLOR=green>" .
                  "<INPUT TYPE='TEXT' NAME='eye_rate$i' READONLY SIZE='1' " .
                  "VALUE='$ByEye[$i]'  STYLE='font-family:courier;font-size:7pt;color:blue; " .
                  "background-color:gold'></TD>" . 

                  "<TD><SELECT NAME='eye$i' STYLE='font-family:courier; font-size:7pt;" .
                  " color:blue; background-color:gold' " .
                  #"OnChange=\"UpdateEye(eye_rate$i,eye$i.value,$RNAMotifMatches[$i]{'UID'})\"  "  .
                  "><OPTION VALUE='-'>Select Score" .
                  "<OPTION VALUE='0'>0: Excellent" .
                  "<OPTION VALUE='1'>1: Good" .
                  "<OPTION VALUE='2'>2: OK" .
                  "<OPTION VALUE='3'>3: Average" .
                  "<OPTION VALUE='4'>4: Poor" .
                  "<OPTION VALUE='5'>5: Bad" .
                  "</SELECT>" .
                  "</TD></TR>\n</TABLE>\n";
     ### END UPDATE BOX ### 

     ### START LINKS IMAGES ###

     $LinksImages = "\n<TABLE CELLPADDING=1 CELLSPACING=1 " .
                    "STYLE='border:1px;border-style:solid;background-color:white'>" . 
                    "\n<TR><TD>" .
                    "\n<A HREF='http://$IPs/cgi-rna/GeneMapper.pl?" . 
                    "Table=$Table&evalue=0.01&QueryText=" .
                    "$RNAMotifMatches[$i]{'CG'}'>" .
                    "<IMG SRC='http:///$IP/Images/map.gif' BORDER=0></A>" .

                    "\n<A HREF='http://$IPs/cgi-rna/ConservationAssessment.pl?" . 
                    "TB=$Table&Sequence=$RNAMotifMatches[$i]{'Sequence'}'>" .
                    "<IMG SRC='http://$IP/Images/miniscales3.gif' BORDER=0></A>";

     if( $RNABinding{$RNAMotifMatches[$i]{'CG'}} eq 1 )
       {
         $LinksImages .= "\n<A HREF='http://$IPs/cgi-rna/RNABindingDB.pl?Phase=Individual" .
                         "&CG=$RNAMotifMatches[$i]{'CG'}&Table=$Table'>" .
                         "<IMG SRC='http://$IP/Images/magnet_small.gif' BORDER=0></A>";
       }
     else
        {
          $LinksImages .= "\n<IMG SRC='http://$IP/Images/blank_small.gif' BORDER=0>";
        }

     if( $Localizing{$RNAMotifMatches[$i]{'CG'}} eq 1 )
        {
          $LinksImages .= "\n<A HREF='http://$IPs/cgi-rna/LocalizationAssessment.pl?Phase=" .
                          "ViewAll#$RNAMotifMatches[$i]{'CG'}&Table=$Table'>" .
                          "<IMG SRC='http://$IP/Images/train.jpg' BORDER=0>";
        }
     else
        {
          $LinksImages .= "\n<IMG SRC='http://$IP/Images/blank_small.gif' BORDER=0>";
        }
     $LinksImages .= "</TD></TR>" .
                     "\n<TR><TD ALIGN=center><INPUT TYPE=button NAME='BTN_Desc_$i' " . 
                     "VALUE='View RNAMotif' " .
                     "onClick='OpenStructureWindowNew(\"http://$IP/RNAMotifs/$In_RNAMotifDesc\");" .
                     "return false' " .
                     "STYLE='font-size:6pt;color:blue;background-color:gold'></TD>" . 
                     "\n</TR></TABLE>";

     ### END LINKS IMAGES ###

     ### START TOP PANEL ###

     $O_Top = "\n<TABLE WIDTH=100% " .
              "STYLE='border:1px;border-style:solid;border-color:white;background-color:white'>\n" .

              "\n<TR>" . 
              "\n<TD WIDTH=50 BGCOLOR=silver><FONT FACE=courier SIZE=1>UID:" . 
              "<FONT COLOR=red>$RNAMotifMatches[$i]{'UID'}</FONT><BR>" . 
              "Rank:<FONT COLOR=red>$Rank</FONT></TD>" . 

              "\n<TD BGCOLOR=silver><FONT FACE=sans,arial SIZE=1>Gene:$RNAMotifMatches[$i]{'CG'} " . 
              "&nbsp Chr:$RNAMotifMatches[$i]{'Chromosome'}" .
              "<BR><INPUT TYPE='TEXT' NAME='Annot$i' READONLY SIZE='40' " .
              "VALUE='$RNAMotifMatches[$i]{'Annotation'}' " . 
              "STYLE='font-family:courier; font-size:7pt; color:blue; " . 
              "background-color:gold'>" . 
              "<BR>RNAdistance Scores grk:[$DistanceScores[0]] ifactor:[$DistanceScores[2]]" .
              "</TD>" .

              "\n<TD BGCOLOR=silver ALIGN=center >$UpdateBox</TD>" .
              "\n<TD BGCOLOR=silver VALIGN=top ALIGN=center>$LinksImages" . 
              "</TD>\n</TR>\n</TABLE>\n";

     ### END TOP PANEL ###

     ### START MID PANEL ###
     $O_Mid =  "<TABLE WIDTH=100% " . 
               "STYLE='border:1px;border-style:solid;border-color:white;background-color:white'>";

     
     $O_Mid .= "<TR>" . 
               "<TD BGCOLOR=silver><FONT FACE='sans,arial' SIZE=1>" . 
               "Filter 1: Low Complexity &nbsp" . 
               "<FONT COLOR=blue>$RNAMotifMatches[$i]{'Complexity'}</FONT></FONT></TD>" . 
               "<TD BGCOLOR=silver><FONT FACE='sans,arial' SIZE=1>" . 
               "Filter 2: Pseudoobscura &nbsp" . 
               "<FONT COLOR=blue>xxx</FONT></TD>" . 
               "<TD BGCOLOR=silver><FONT FACE='sans,arial' SIZE=1>" . 
               "Filter 3: Location &nbsp" . 
               "<FONT COLOR=blue>xxx</FONT></TD>" . 
               "</TR>" .
               "</TABLE>";

     ### END MID PANEL ###

     ### START SEQUENCES PANEL ###     

     if( $RNAMotifMatches[$i]{'RNAMotif_SS'} =~ m/xxx/)
       {
         @RNAMotif_SSs = split(/xxx/,$RNAMotifMatches[$i]{'RNAMotif_SS'});
       }
     else
       {
         $RNAMotif_SSs[0] = $RNAMotifMatches[$i]{'RNAMotif_SS'};
       }

     $RNAMotifDropDown = "<SELECT NAME='RNAMotifs$i' STYLE='font-family:courier;font-size:7pt;" .
                         " color:blue;background-color:gold' " .
                         ">";
     for($j=0; $j<=$#RNAMotif_SSs; $j++)
        {
          $RNAMotifDropDown .= "<OPTION VALUE='$RNAMotif_SSs[$j]'>$RNAMotif_SSs[$j]";
        }
     $RNAMotifDropDown .= "</SELECT>";
     $RNAMotifCount = $#RNAMotif_SSs + 1;

     undef(@RNAMotif_SSs);

     $SeqPanel = "\n<TABLE BGCOLOR=silver STYLE='padding-left:5px;padding-right:5px' CELLSPACING=0 WIDTH=695 " . 
                 "STYLE='border:1px;border-style:solid;border-color:white;background-color:white'>" .
                 "\n<TR><TD WIDTH=55><FONT FACE=courier SIZE=1>Sequence</TD>" . 
                 "\n<TD><FONT FACE=courier SIZE=1>&nbsp</FONT>" . 
                 "<FONT FACE=courier SIZE=1 COLOR=green>$RNAMotifMatches[$i]{'Sequence'}</TD>" .
                 "\n<TD WIDTH=100><FONT FACE=courier SIZE=1 COLOR=red>xxx</FONT></TD>" .
                 "<TD><FONT FACE=courier SIZE=1 COLOR=red>xxx</FONT></TD></TR>\n" .

                 "\n<TR><TD><FONT FACE=courier SIZE=1>Mfold</TD>" . 
                 "\n<TD><FONT FACE=courier SIZE=1>&nbsp<FONT COLOR=blue>$MFOLD</TD>" .
                 "\n<TD><FONT FACE='sans.arial' SIZE=1>&Delta;G = $MEnergy</TD>" .
                 "\n<TD ALIGN=left><INPUT TYPE=button NAME='BTN_MFOLD_$i' VALUE='View MFOLD' " .
                 "onClick='OpenStructureWindowNew(\"/cgi-rna/Version_3/OnTheFlyStructureGIF_v3.pl?" . 
                 "Phase=MFOLD&Desc=$Table&UID=$RNAMotifMatches[$i]{'UID'}" . 
                 "&Seq=$RNAMotifMatches[$i]{'Sequence'}&Size=50\");return false' " .
                 "STYLE='font-size:6pt;color:blue;background-color:gold'></TD></TR>\n" .

                 "\n<TR><TD><FONT FACE=courier SIZE=1>RNAlfold</TD>" . 
                 "\n<TD><FONT FACE=courier SIZE=1>&nbsp<FONT COLOR=blue>$LFOLD</TD>" .
                 "\n<TD><FONT FACE='sans.arial' SIZE=1>&Delta;G = $LEnergy</TD>" .
                 "\n<TD><INPUT TYPE=button NAME='BTN_LFOLD_$i' VALUE='View LFOLD' " .
                 "onClick='OpenStructureWindowNew(\"/cgi-rna/Version_3/OnTheFlyStructureGIF_v3.pl?" .
                 "Phase=LFOLD&Desc=$Table&UID=$RNAMotifMatches[$i]{'UID'}" .
                 "&Seq=$RNAMotifMatches[$i]{'Sequence'}&Size=50\");return false' " .
                 "STYLE='font-size:6pt;color:blue;background-color:gold'></TD></TR>\n" .

                 "\n<TR><TD VALIGN=top><FONT FACE=courier SIZE=1>RNAMotif</TD>" . #<BR>" . 
                 # "<FONT FACE='sans,arial'>S=$RNAMotifMatches[$i]{'RNAMotif_Score'}</FONT></TD>" . 

                 "\n<TD ALIGN=left>$RNAMotifDropDown</TD>" .
                 "\n<TD><FONT FACE='sans.arial' SIZE=1>No. Variants = $RNAMotifCount</TD>" .

                 "\n<TD VALIGN=middle><INPUT TYPE=button NAME='BTN_RNAMotif_$i' " . 
                 "VALUE='View RNAMotif Fold' " .
                 "onClick='OpenStructureWindowNewTest(\"/cgi-rna/Version_3/" . 
                 "OnTheFlyStructureGIF_v3.pl?" .
                 "Phase=RNAMOTIF&Table=R&Desc=$Table&UID=$RNAMotifMatches[$i]{'UID'}&Str=\"" .
                 ",Myform.RNAMotifs$i.options[Myform.RNAMotifs$i.selectedIndex].value," . 


                 "\"&Seq=$RNAMotifMatches[$i]{'Sequence'}&Size=50\");return false' " .
                 "STYLE='font-size:6pt;color:blue;background-color:gold'></TD>" .
                 "</TR></TABLE>\n";

#document.Myform.RNAMotifs$i.options[document.Myform.RNAMotifs$i.selectedIndex].value

     ### END SEQUENCES PANEL ###

     $Output .= "\n<TABLE WIDTH=700 " . 
                "STYLE='border:1px;border-style:solid;background-color:white'>" . 

                "\n<TR><TD>$O_Top</TD></TR>" .
                "\n<TR><TD>$O_Mid</TD></TR>" .
                "\n<TR><TD>$SeqPanel</TD>" . 
                "</TR>\n</TABLE>\n<P>";

     undef($O_Top);
     undef($O_Left);
     undef($UpdateBox);
   }

$dbh->disconnect;
undef(@RNAMotifMatches);
return $Output;
}

#------------------------------------------------------------------------------
sub PrintOutput_Pass  {
                                                                                        
my $Output = $_[0];
                                                                                        
print "Content-type: text/html\n\n";
                                                                                        
print <<END_OF_PAGE;
                                                                                        
<HTML>
<HEAD>
<TITLE>RNAMotif Search Results: Russell S. Hamilton</TITLE>
                                                                                        
<SCRIPT LANGUAGE="JavaScript">
<!--hide
                                                                                        
                                                                                        
function OpenStructureWindow(link)
   {
      StructureWindow = window.open(link,"StructureWindow", "toolbar=no, location=no, menubar=no, scrollbars=yes, width=620, height=800");
                                                                                        
      StructureWindow.blur();
   }
                                                                                        
function OpenStructureWindowNew(link)
   {
      StructureWindowNew = window.open(link,"StructureWindowNew", "toolbar=no, location=no, menubar=no, scrollbars=yes, width=460, height=700");
                                                                                        
                                                                                        
      StructureWindowNew.blur();
   }

function OpenStructureWindowNewTest(link, link2, link3)
   {
      StructureWindowNew = window.open(link+link2+link3,"StructureWindowNew", "toolbar=no, location=no, menubar=no, scrollbars=yes, width=460, height=700");

      StructureWindowNew.blur();
   }

function Test(ID)
   {
     var Value;

     Value = ID + "dave";

     return Value;
   }

function UpdateCandidate(can,id)
   {
      var candidate;
      var link;
                                                                                        
      if(can.checked == true)
        {
          candidate = 1;
        }
      else
        {
          candidate = 0;
        }
                                                                                        
      link = "http://$IPs/cgi-rna/UpdateCandidates.pl?CANDIDATE="
             + candidate + "&ID=" + id + "&Table=$Table";
                                                                                        
      CandidateWindow = window.open(link, "CandidateWindow", "toolbar=no, location=no, menubar=no, scrollbars=yes, width=100, height=100");
      CandidateWindow.close();
   }
                                                                                        
function UpdateComments(com,id)
   {
      var comments;
      var link;
                                                                                        
      comments  = com.value;
                                                                                        
      link = "http://$IPs/cgi-rna/UpdateComments.pl?COMMENTS="
             + comments + "&ID=" + id + "&Table=$Table";
                                                                                        
                                                                                        
      CommentsWindow = window.open(link, "CommentsWindow", "toolbar=no, location=no, menubar=no, scrollbars=yes, width=100, height=100");
      CommentsWindow.close();
   }
                                                                                        
function UpdateEye(object, eye,id)
   {
     var link;
                                                                                        
     link = "http://$IPs/cgi-rna/UpdateEye.pl?EYE=" + eye + "&ID=" + id + "&Table=$Table";
                                                                                        
     object.value = eye;
                                                                                        
     EyeWindow = window.open(link, "EyeWindow", "toolbar=no, location=no, menubar=no, scrollbars=yes, width=100, height=100");
     EyeWindow.close();
   }
                                                                                        
                                                                                        
//-->

</SCRIPT>
                                                                                        
</HEAD>
<BODY BGCOLOR=#9999FF><FONT FACE="helvetica" SIZE=2>
@Header <P>
<FONT FACE=courier COLOR=white SIZE=5>RNAMotif Results</FONT><BR>
                                                                                        
<FORM NAME="Myform" method="post" action="">
$Output
                                                                                        
</FORM>
<P>@Tail
</BODY>
</HTML>
END_OF_PAGE
}

#------------------------------------------------------------------------------
# FIN
#------------------------------------------------------------------------------
