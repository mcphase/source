#!/usr/bin/perl
# rhoso1ion - program to calculate simple impurity scattering resistivity from crystal field

BEGIN{@ARGV=map{glob($_)}@ARGV}

use Getopt::Long;
use Math::Complex;

$SMALL=1e-5;
#$debug = true;

# -------------------------------------------------------------------------------------- #
# Subfunctions: 
#   mmult / range / veclen / matdim  -  Matrix multiplication from Perl Cookbook
#   zeros - Creates a zero matrix
#   usage - Prints help screen
# -------------------------------------------------------------------------------------- #

# The following subroutines were shamelessly stolen from the Perl Cookbook
# by O'Reily (Recipe 2.14), and used for the matrix multiplications.

sub mmult {
    my ($m1,$m2) = @_;
    my ($m1rows,$m1cols) = matdim($m1);
    my ($m2rows,$m2cols) = matdim($m2);
    unless ($m1cols == $m2rows) {       # raise exception
        die "IndexError: matrices don't match: $m1cols != $m2rows";
    }
    my $result = [];
    for $i (range($m1rows)) { for $j (range($m2cols)) { $result->[$i][$j]=0; } }
    my ($i, $j, $k);
    for $i (range($m1rows)) {
        for $j (range($m2cols)) {
            for $k (range($m1cols)) {
                $result->[$i][$j] += $m1->[$i][$k] * $m2->[$k][$j];
            }
        }
    }
    return $result;
}

sub range { 0 .. ($_[0] - 1) }

sub veclen {
    my $ary_ref = $_[0];
    my $type = ref $ary_ref;
    if ($type ne "ARRAY") { die "$type is bad array ref for $ary_ref" }
    return scalar(@$ary_ref);
}

sub matdim {
    my $matrix = $_[0];
    my $rows = veclen($matrix);
    my $cols = veclen($matrix->[0]);
    return ($rows, $cols);
}

# End of stolen code!

sub zeros {
    my ($rows,$cols) = @_;
    my $ii,$jj;
    my @result;
    for $ii(0..($rows-1)) { my @blank; for $jj(0..($cols-1)) { push(@blank,0); } push(@result,\@blank); }
    return \@result;
}

sub usage() {
  print STDERR << "EOF";
    $0:
    Calculates the impurity scattering resistivity due to the exchange
    interaction of a conduction s-electron with a localised f-electron
    in a crystal field, from the output of the program so1ion.
    The theory used in this calculation is given in the source code of
    this program, and is based on equations 6, 7, 8 of Rao and Wallace
    PRB vol 2, pp 4613 (1970) and 8 and 9 of Dekker, J Appl Phys vol36
    pp 906 (1965). 

    usage: $0 Tmin Tmax deltaT
    usage: $0 col1 col2 datafile -r rho0

    Same syntax as for cpso1ion. The second option will output a sta
    value, but you should use the option -r or --rho0 to specify the
    factor proportional to the s-conduction electron exchange.

    usage: $0 [-h] [--help]
              [-t <T_list>] [--temperature <T_list>]
                  (default is 1 to 300K in 1K steps)
              [-f <T_file>] [--file <T_file>]
              [-r <rho_0>] [--rho0 <rho_0>]
              [-i CEF_matrix_file] [--input CEF_matrix_file]
                  (default is results/levels.cef)

    The T_list can be given as a list of comma separated numbers, 
    T1,T2,T3,... or a set of 3 colon separated numbers: start:step:end 
    Alternatively the temperature at which to calculate can be given
    in a file. If this is followed by a column number, that column will
    be assumed to contain the temperatures to calculate at (otherwise
    the first column will be used). If there is another column number,
    this column will be overwritten with the calculated resistivity. 
    Otherwise the program outputs the results to the screen (standard
    output). 

    Note that if the options are given, parameters in them (e.g. Tmin
    Tmax deltaT) take precedent over non-option parameters.

    The calculated resistivity is multiplied by rho_0 which is equal
    to (3pi N m/hbar e^2 E_f) * G^2 * (g-1)^2, where N is the number
    of f-electrons, m the electron mass, G the s-f exchange and g the
    Lande g-factor.

    E.g.  $0 1 4 expdat.res
                 will take temperature from file expdat.res column 1
                 and calculate the difference between calculated
                 and experimental resistivity in column 4.
          $0 -f expdat.res 1 4
                 will take temperature from file expdat.res column 1
                 and write calculated resistivity to column 4
          $0 1 300 0.5 > calc.res
          $0 -T 1:0.5:300 > calc.res
                 will calculate the resistivity from 1 to 300K in 0.5K
                 steps and output to file calc.res
EOF
  exit;

}

sub printheader {
  my $filehandle;
  if($#_==-1) { $filehandle=STDOUT; } else { $filehandle=$_[0]; }
  print { $filehandle } "#! J=",$J," Eigenvalues: ",join("\t",@E),"\n";
  print { $filehandle } "#! Eigenvectors: "; 
  for($ii=0; $ii<=$J2; $ii++) { for($jj=0; $jj<=$J2; $jj++) {
    print { $filehandle } $V->[$ii][$jj]*1,"\t"; } print { $filehandle } "\n#!\t\t",; } 
  print { $filehandle } "\n";
}


# see http://aplawrence.com/Unix/perlgetops.html for details of GetOptions

usage() if $#ARGV<0;

# -------------------------------------------------------------------------------------- #
# Options and declarations
# -------------------------------------------------------------------------------------- #
GetOptions("help"=>\$helpflag,
           "rho0=f"=>\$rho0,
           "temperature=s"=>\$tlist,
	   "file=s"=>\$tfile,
	   "input=s"=>\$input);

usage() if $helpflag;

$k_B = 1.3806505e-23; $Q_e = 1.60217653e-22;
$i=Math::Complex->make(0,1);
$Tmin=1; $Tmax=300; $deltaT=1; $col1=1;   # Defaults

# -------------------------------------------------------------------------------------- #
# Parse input
# -------------------------------------------------------------------------------------- #

# Parses the temperature / file inputs
if ($#ARGV==2) {     # Three input arguments besides flags.
  if (-e $ARGV[2]) { $col1=$ARGV[0]; $col2=$ARGV[1]; $datfile=$ARGV[2]; }
  else { $Tmin=$ARGV[0]; $Tmax=$ARGV[1]; $deltaT=$ARGV[2]; }
}
elsif ($#ARGV==1) { $col1=$ARGV[0]; $col2=$ARGV[1]; }

if ($tlist) { 
     if ($tlist =~/,/) { @Temp = split(/,/,$tlist); }
  elsif ($tlist =~/:/) { @_ = split(/:/,$tlist); $Tmin=$_[0]; $deltaT=$_[1]; $Tmax=$_[2]; }
}

if ($tfile) { if (-e $tfile) { open(INFILE,$tfile); @Temp=(); $fflag=1; }
  else { warn "$0: File $tfile does not exist. Ignoring."; } }
elsif ($datfile) { if (-e $datfile) { open(INFILE,$datfile); @Temp=(); $fflag=1; $dflag=1; }
  else { warn "$0: File $datfile does not exist. Ignoring."; } }

if($fflag) {
  while(<INFILE>) { 
    push(@INLines,$_);
    if ($_!~/^\s*#/) { 
      split; push(@Temp,$_[$col1-1]); if($col2 && $dflag) { push(@resexp,$_[$col2-1]); }
    } 
  }
  close(INFILE);
}

if (!@Temp) { for ($iT=$Tmin; $iT<=$Tmax; $iT+=$deltaT) { push(@Temp,$iT); } }

# Reads in CF energies and wavefunctions
if (!$input) { $ceflevels = "results/levels.cef"; } else { $ceflevels=$input; }
open (CEFLEV, $ceflevels) or die "$0: cannot open $ceflevels for input CF eigenvalues/vectors.";
$rpart=0; $ipart=0; $einit=0;
while(<CEFLEV>) {
  if ($_=~/^\s*J=/ || $_=~/^#!J=/) { @Jline=split; $J=$Jline[1]; next; } #$V=zeros(2*$J+1,2*$J+1);
  if ($_=~/Eigenvalues/) { $_ =~ tr/0-9\ \.\-//cd; @ereal=split; $JreadR=($#ereal)/2; $einit=1; }
  if ($_=~/^\s*#/) { 
    if ($_=~/Real/)      { $rpart=1; $ipart=0; $ct=$einit; }
    if ($_=~/Imaginary/) { $rpart=0; $ipart=1; $ct=$einit; }
  }
  else {
    if ($rpart) {
      if (!$J && $ct==0) { warn("$0: Line 'J=' not found before '#Real part' in $ceflevels file."); } 
      if ($ct==0) {
        @ereal=split; $JreadR=($#ereal)/2; $ct++;
        if ($J && $JreadR!=$J) { warn("Number of real eigenvalues != 2J+1"); }
      }
      else { push(@vreal,split); }
    }
    if ($ipart) {
      if (!$J && $ct==0) { warn("$0: Line 'J=' not found before '#Imaginary part' in $ceflevels file."); } 
      if ($ct==0) {
        @eimag=split; $JreadI=($#eimag)/2; $ct++;
        if ($J && $JreadI!=$J) { warn("Number of imaginary eigenvalues != 2J+1"); }
      }
      else { push(@vimag,split); }
    }
  }   
}
close(CEFLEV);

if ($einit==0 && $JreadR!=$JreadI) { 
  die "$0: Real and Imaginary eigenvector matrix do not have the same dimensions."; }
$J=$JreadR; $J2 = 2*$J; $Jsq=$J*($J+1);

# -------------------------------------------------------------------------------------- #
# Calculates the angular momentum operator matrices.
# -------------------------------------------------------------------------------------- #
$szJz = zeros(2*($J2+1),2*($J2+1));  # Matrix sz.Jz between s conduction and f electrons
$spJm = zeros(2*($J2+1),2*($J2+1));  # Matrix s+.J- between s conduction and f electrons
$smJp = zeros(2*($J2+1),2*($J2+1));  # Matrix s-.J+ between s conduction and f electrons

for $ii(0..($J2)) { push(@mj,$ii-$JreadR); $szJz->[$ii][$ii]=$mj[-1]/2; $szJz->[$ii+$J2+1][$ii+$J2+1]=-$mj[-1]/2; }
push(@mj,@mj);
for $ii(0..(2*$J2+1)) { for $jj(0..(2*$J2+1)) {
#    if(  ($mj[$jj]==$mj[$ii]-1) ) { $spJm->[$ii][$jj] = sqrt($Jsq-$mj[$ii]*($mj[$ii]-1)); }
# elsif(  ($mj[$jj]==$mj[$ii]+1) ) { $smJp->[$ii][$jj] = sqrt($Jsq-$mj[$ii]*($mj[$ii]+1)); }
     if( ($ii-$jj)>($J2) && ($mj[$jj]==$mj[$ii]-1) ) { $spJm->[$ii][$jj] = sqrt($Jsq-$mj[$ii]*($mj[$ii]-1)); }
  elsif( ($jj-$ii)>($J2) && ($mj[$jj]==$mj[$ii]+1) ) { $smJp->[$ii][$jj] = sqrt($Jsq-$mj[$ii]*($mj[$ii]+1)); }
} }

$sJmat=[]; for $ii(0..(2*$J2+1)) { for $jj(0..(2*$J2+1)) {
  $sJmat->[$ii][$jj] = $szJz->[$ii][$jj] + ($spJm->[$ii][$jj] + $smJp->[$ii][$jj])/2; 
} }

#ms = [ones(1,sz(1))/2 -ones(1,sz(1))/2]; mj = [-J:J -J:J]; HszJz = diag(ms).*diag(mj); 
#%HspJm = zeros(2*sz(1)); HsmJp = zeros(2*sz(1));
#%for i=1:(2*sz(1))
#%  for j=1:(2*sz(2))
#%        if(ms(j)==(ms(i)+1) && mj(j)==(mj(i)-1)); HspJm(i,j) = sqrt(J*(J+1)-mj(i)*(mj(i)-1)); 
#%    elseif(ms(j)==(ms(i)-1) && mj(j)==(mj(i)+1)); HsmJp(i,j) = sqrt(J*(J+1)-mj(i)*(mj(i)+1)); end
#%  end
#%end
#HsmJp = [zeros(sz) Jmat(:,:,5).*2; zeros(sz) zeros(sz)];
#HspJm = [zeros(sz) zeros(sz); Jmat(:,:,4).*2 zeros(sz)];
#
#sJmat = HszJz + (HspJm+HsmJp)/2;

if ($debug) {
  print "s.J matrix:\n";
  for $ii(0..(2*$J2+1)) { for $jj(0..(2*$J2+1)) {
    printf("%0.5g\t",$sJmat->[$ii][$jj]); } print "\n"; }
}

# -------------------------------------------------------------------------------------- #
# Parses energy/wavefunction and calculates the Bose and Fermi functions.
# -------------------------------------------------------------------------------------- #

# Checks that the eigenvector matrix is real
$isreal=1; foreach (@vimag) { if (abs($_)>$SMALL) { $isreal=0; break; } }
foreach (@eimag) { if (abs($_)>$SMALL) { warn("Imaginary eigenvalues (energies). Ignoring."); break; } }

# Populates the Energy and Wavefunction (Eigenvalue/vector) arrays.
$V=zeros($J2+1,$J2+1); @E=@ereal; $ix=0;
if($isreal) {
  for ($ii=0; $ii<($J2+1); $ii++) { for ($jj=0; $jj<($J2+1); $jj++) {
    $V->[$ii][$jj] = $vreal[$ix++]; } }
}
else {
  for ($ii=0; $ii<($J2+1); $ii++) { for ($jj=0; $jj<($J2+1); $jj++) {
    $V->[$ii][$jj] = Math::Complex->make($vreal[$ix],$vimag[$ix++]); } }
}
$Vi=zeros(2*$J2+2,2*$J2+2); $tVi=zeros(2*$J2+2,2*$J2+2);
for $ii(0..(2*$J2+1)) { for $jj(0..(2*$J2+1)) {
  $Vi->[$ii][$jj] = $V->[$ii%($J2+1)][$jj%($J2+1)]; $tVi->[$jj][$ii] = $Vi->[$ii][$jj]; } }

# Calculates the matrix elements |<V|s.J|V>|^2
$sJVi = mmult($sJmat,$Vi); $matel = mmult($tVi,$sJVi);
for $ii(0..(2*$J2+1)) { for $jj(0..(2*$J2+1)) { 
  if($isreal) { $matel->[$ii][$jj]*=$matel->[$ii][$jj]; } else { $matel->[$ii][$jj]*=~$matel->[$ii][$jj]; }
} }

if ($debug) {
  print "Eigenvalues:\n",join("\t",@E),"\n";
  print "Eigenvectors:\n"; 
  for($ii=0; $ii<=$J2; $ii++) { for($jj=0; $jj<=$J2; $jj++) {
    print $V->[$ii][$jj]*1; print "\t"; } print "\n"; }
  print "Extended Eigenvectors:\n"; 
  for $ii(0..(2*$J2+1)) { for $jj(0..(2*$J2+1)) {
    print $Vi->[$ii][$jj]*1; print "\t"; } print "\n"; }
  print "V'*sJ*V:\n"; 
  for $ii(0..(2*$J2+1)) { for $jj(0..(2*$J2+1)) {
    printf("%0.4g\t",$matel->[$ii][$jj]); } print "\n"; }
}

# Converts energy from meV to K to save computation later
for ($ii=0; $ii<=$J2; $ii++) { push(@Ek,$E[$ii]*$Q_e/$k_B); }

# Calculates pi and fij
foreach $T (@Temp) {
  $Z=0; for ($ii=0; $ii<=$J2; $ii++) { $Z += exp(-$Ek[$ii]/$T); }
  @p=(); @f=();
  for ($ii=0; $ii<=$J2; $ii++) {
    push(@p,exp(-$Ek[$ii]/$T)/$Z);
    for ($jj=0; $jj<=$J2; $jj++) {
      push(@f,2/(1+exp(-($Ek[$ii]-$Ek[$jj])/$T)));
  } }

  if ($debug && !$firstloop) {
    print "Z = $Z\n"; $firstloop=1;
    print "p_i = "; for ($ii=0; $ii<=$J2; $ii++) { printf("%0.3g\t",$p[$ii]); } print "\n";
    print "f_ij=\t";for ($ii=0; $ii<=$J2; $ii++) { for ($jj=0; $jj<=$J2; $jj++) {
      printf("%0.3g\t",$f[$jj+$ii*($J2+1)]); } print "\n"; if($ii<$J2) { print "\t"; } }
  }

  $res=0;
  for $ii(0..(2*$J2+1)) { for $jj(0..(2*$J2+1)) {
    $xi=$ii%($J2+1); $res += $matel->[$ii][$jj] * $p[$xi] * $f[$jj%($J2+1)+$xi*($J2+1)]; } }

  $res /= $J*($J+1)*(4/3); 
  if($rho0) { push(@reslist,$rho0*$res); } else { push(@reslist,$res); }
# print "$T\t$res\n";

#Vi = [V V; V V]; res = zeros(size(T)); p = [p p]; f = [f f; f f];
#for i=1:length(Vi)
#  for j=1:length(Vi)
#    res = res + ( Vi(:,j)' * sJmat * Vi(:,i) ).^2.*p{i}.*f{i,j};
#  end
#res = res ./ (J*(J+1));
#res = res ./ (4/3);      % Unknown correction...

}


# -------------------------------------------------------------------------------------- #
# Outputs
# -------------------------------------------------------------------------------------- #

# User specified a temperature datafile and two columns. Replace column2 with calculation.
if ($tfile && $fflag && $col2) {
  open($OUTFILE,">rhoso1ion.out"); $ix=0;
  printheader($OUTFILE);
  foreach(@INLines) {
    if ($_=~/^\s*#/) { print $OUTFILE $_; }
    else {
      split; if($#_<($col2-1)) { for $ic ($#_..($col2-2)) { push(@_,0); } }
      $_[$col2-1] = $reslist[$ix++];
      print $OUTFILE join("\t",@_),"\n";
    }
  }
  close($OUTFILE); 
  unless(rename("rhoso1ion.out",$tfile)) {
    open ($OUTFILE, ">$tfile") or die "$0: cannot to write to $tfile. Output left in rhoso1ion.out";
    open($INFILE,"rhoso1ion.out");
    while(<$INFILE>) { print $OUTFILE $_; }
    close($INFILE); close($OUTFILE); unlink("rhoso1ion.out");
  }
  exit(0);
}

# User specified a datafile with resistivity data, and wants sta outputted.
if ($fflag && $dflag) {
  printheader();
  $sta=0; print "# T(K)  Res_Calc  Res_exp\n";
  for (0..$#reslist) {
    $sta+=($reslist[$_]-$resexp[$_])**2;
    print $Temp[$_]," ",$reslist[$_]," ",$resexp[$_],"\n";
  }
  $sta /= ($#reslist+1);
  print "sta=$sta\n";
  exit(0);
}

printheader();
print "# T(K)  Res_Calc\n";
for (0..$#reslist) { print $Temp[$_]," ",$reslist[$_],"\n"; }

# -------------------------------------------------------------------------------------- #
# Original Matlab Code
# -------------------------------------------------------------------------------------- #
#function res = cfres(T,Hcf)
#
#k_B = 1.3806505e-23; Q_e = 1.60217653e-22;
#
#sz = size(Hcf);
#if sz(1)~=sz(2); error('Hcf is not square!'); end
#
#J = (sz(1)-1)/2;
#Jmat = mag_op_j(J);
#
#[V,E] = eig(Hcf); E = diag(E)-min(min(E));
#
#Z = zeros(size(T)); for i=1:length(E); Z = Z + exp(-E(i).*Q_e./(k_B.*T)); end
#for i=1:length(E)
#  p{i} = exp(-E(i).*Q_e./(k_B.*T))./Z; 
#  for j=1:length(E)
#   %f{i,j} = 2./(1 + exp(-abs(E(i)-E(j)).*Q_e./(k_B.*T)));
#    f{i,j} = 2./(1 + exp(-(E(i)-E(j)).*Q_e./(k_B.*T)));
#  end
#end
#
#% rho_s(T)/rho_0 = sum_{ms,ms',i,i'} <ms',i'| szJz + (1/2)(s+J- + s-J+) | ms,i>^2 p_i f_ii'  /  J(J+1)  [Rao and Wallace]
#% <ms,mj| szJz | ms,mj> = msmj      |     <ms+/-1,mj-/+1| s+/-J-/+ | ms,mj> = sqrt(J(J+1)-mj(mj-/+1))   [Dekker]
#
#ms = [ones(1,sz(1))/2 -ones(1,sz(1))/2]; mj = [-J:J -J:J]; HszJz = diag(ms).*diag(mj); 
#%HspJm = zeros(2*sz(1)); HsmJp = zeros(2*sz(1));
#%for i=1:(2*sz(1))
#%  for j=1:(2*sz(2))
#%        if(ms(j)==(ms(i)+1) && mj(j)==(mj(i)-1)); HspJm(i,j) = sqrt(J*(J+1)-mj(i)*(mj(i)-1)); 
#%    elseif(ms(j)==(ms(i)-1) && mj(j)==(mj(i)+1)); HsmJp(i,j) = sqrt(J*(J+1)-mj(i)*(mj(i)+1)); end
#%  end
#%end
#HsmJp = [zeros(sz) Jmat(:,:,5).*2; zeros(sz) zeros(sz)];
#HspJm = [zeros(sz) zeros(sz); Jmat(:,:,4).*2 zeros(sz)];
#
#sJmat = HszJz + (HspJm+HsmJp)/2;
#Vi = [V V; V V]; res = zeros(size(T)); p = [p p]; f = [f f; f f];
#for i=1:length(Vi)
#  for j=1:length(Vi)
#    res = res + ( Vi(:,j)' * sJmat * Vi(:,i) ).^2.*p{i}.*f{i,j};
#  end
#end
#
#res = res ./ (J*(J+1));
#
#res = res ./ (4/3);      % Unknown correction...
#
#%me = zeros(size(T));
#%for i=1:length(Vi)
#%  me = me + (Vi(:,j)' * sJmat * Vi(:,i)).*p{i};
#%end
#%plot(1:300,me./Z);
