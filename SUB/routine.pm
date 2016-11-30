#package routine;
use strict;

sub FBPATHS {
  my $location = "/home/modupe/public_html/FBTranscriptAtlas/";
  our $chickenpath = $location."chickenvariants";
  our $mousepath = $location."mousevariants";
  our $alligatorpath = $location."alligatorvariants";
  return ($chickenpath, $mousepath, $alligatorpath);
  
}
sub GETINFO {
  my $tisue = uc $_[0];
  my %GENES = %{$_[1]};
  print "<xx>$tisue<xx>";
    print "<table class=\"gened\">
            <tr>
              <th class=\"gened\">Line</th>
              <th class=\"gened\">Maximum Fpkm</th>
              <th class=\"gened\">Average Fpkm</th>
              <th class=\"gened\">Minimum Fpkm</th>
            </tr>\n";
      foreach my $a (keys %GENES){
        print "<tr>
              <th class=\"geneds\" colspan=100%>$a</th></tr>\n";
        foreach my $b (sort keys % {$GENES{$a} }){
          my @all = split('\|', $GENES{$a}{$b}, 3);
          print "<tr><td class=\"gened\"><b>$b</b></td>
                <td class=\"gened\">$all[0]</td>
                <td class=\"gened\">$all[1]</td>
                <td class=\"gened\">$all[2]</td></tr>\n";
        }    
      }
      print "</table>\n";
      undef %GENES;
}
1;
