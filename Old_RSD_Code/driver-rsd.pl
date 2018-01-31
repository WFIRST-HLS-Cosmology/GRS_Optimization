# Driver for RSD code that uses the same input file format as the BAO tool.

$line = <>;
($N, $fsky, $sig) = split ' ', $line;
$A = $fsky * 41253;

$FG = 0;
while ($line=<>) {
  $line =~ s/\n//;
  @L = split ' ', $line;
  $line_in = "$L[0] $L[1] $L[2] $L[3]";
  print "$L[0]\t$L[1]\t$L[2]\t$L[3]\t";
  ($sigf) = split ' ', qx/python rsdmat.py $line_in $A $sig/;
  print "$sigf\n";
  $FG += 1./$sigf**2;
}

print (sprintf "# %7.1f\n", $FG);
