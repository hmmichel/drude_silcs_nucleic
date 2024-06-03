#!/usr/bin/perl

use strict;

# Takes an input PDB file from a SILCS simulation in GROMACS and adds 
# segid to the file. This is for nucleic acids and it makes some hard-
# coded assumptions. Beware!
#
# Example usage:
#   perl write_segid.pl input.pdb output.pdb
#

unless(scalar(@ARGV)==2)
{
    die "Usage: perl $0 input.pdb output.pdb\n";
}

my $input = $ARGV[0];
my $output = $ARGV[1];

my $pdbfmt = "%-6s%5d%5s %-4s%1s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%10s%2s\n";

open(IN, "<$input") || die "Cannot open $input: $!\n";
my @in = <IN>;
close(IN);

# output coordinate file, name specified by user
open(OUT, ">$output") || die "Cannot open $output: $!\n";

# output file to be used by CHARMM - records properties of system
# like the number of each fragment type
open(STROUT, ">system_properties.str") || die "Cannot open system_properties.str: $!\n";
print STROUT "\* system properties for SILCS\n";
print STROUT "\*\n\n";

# counters for atoms and residue renumbering
my $atoms = 0;
my $offset = 0;

# save the DNAA sequence
my $seq   = '';

# counters for keeping track of the number of each molecule/chain length
my $ndnaa = 0;  # DNA chain length - useful for patch loop in CHARMM script
my $nheta = 0;  # number of HETA ions - useful for GQs
my $nsolv = 0;  # number of water molecules
my $nacey = 0;  # acetate (ACEY)
my $nbenx = 0;  # benzene (BENX)
my $ndmee = 0;  # dimethylether (DMEE)
my $nform = 0;  # formamide (FORM)
my $nimid = 0;  # imidazole (IMID)
my $nmamy = 0;  # methylammonium (MAMY)
my $nmeoh = 0;  # methanol (MEOH)
my $nprpx = 0;  # propane (PRPX)

# box information
my $boxx  = 0;
my $boxy  = 0;
my $boxz  = 0;
my $alp   = 0;
my $bet   = 0;
my $gam   = 0;
my $shape = '';

for (my $i=0; $i < scalar(@in); $i++)
{
    # get box size and shape information
    if ($in[$i] =~ /CRYST1/)
    {
        my @tmp = split(" ", $in[$i]);
        $boxx = $tmp[1];
        $boxy = $tmp[2];
        $boxz = $tmp[3];
        $alp  = $tmp[4];
        $bet  = $tmp[5];
        $gam  = $tmp[6];
        
        if ($boxx == $boxy)
        {
            if ($boxy == $boxz)
            {
                $shape = "cubic";
            }
            else
            {
                $shape = "tetragonal";
            }
        }
        else
        {
            if (($alp == 90.0) && ($gam == 90.0))
            {
                if ($bet == 90.0)
                {
                    $shape = "orthorhombic";
                }
                else
                {
                    $shape = "monoclinic";
                }
            }
            else
            {
                $shape = "triclinic";
            }
        }
    }

    # parse atoms
    if ($in[$i] =~ /ATOM/)
    {
        $atoms++;
        # we have an atom entry, so add the segid
        my ($atom, $anum, $name, $resn, $chain, $resnr, $x, $y, $z, $b, $occ) = unpack("A6A5A5A5A1A5A12A8A18A6A6", $in[$i]);

        # detect residue name and decide what to name the corresponding segment
        # nucleic acid defaults to DNAA, fix later if you're doing something else
        # SOL renamed to TIP3 and its atoms will be renamed OW -> OH2, HW1 -> H1, HW2 -> H2
        # any SILCS fragments will be assigned its residue name as segid

        # remove whitespace from residue name that may cause misalignment when writing output PDB
        $resn =~ s/\s//g;
        my $segid = '';
        if ( ($resn =~ /ADE/) || ($resn =~ /CYT/) || ($resn =~ /GUA/) || ($resn =~ /THY/) || ($resn =~ /DA/) || ($resn =~ /DC/) || ($resn =~ /DG/) || ($resn =~ /DT/))
        {
            $segid = "DNAA";
            $resn =~ s/DA/ADE/;
            $resn =~ s/DC/CYT/;
            $resn =~ s/DG/GUA/;
            $resn =~ s/DT/THY/; 
            # keep track of how many residues there are
            if ($name =~ /C5'/)
            {
                $ndnaa++;
                $seq .= "$resn "
            }
        }
        elsif ($resn =~ /IMIA/)
        {
            $resn =~ s/IMIA/IMID/;
            $segid = $resn;
            #print "Found IMIA, using res = $resn and seg = $segid\n";
        }
        elsif ($resn =~ /MEOH/)
        {
            $segid = $resn;
            # rename atoms
            $name =~ s/CB/C1/;
            $name =~ s/OG/O1/;
            $name =~ s/HG1/HO1/;
            $name =~ s/HB1/H1A/;
            $name =~ s/HB2/H1B/;
            $name =~ s/HB3/H1C/;
        }
        elsif ($resn =~ /POT/)
        {
            $segid = "HETA";
            $nheta++;   # one ion per residue in additive FF
        }
        elsif ($resn =~ /SOL/)
        {
            $segid = "SOLV";
            # rename atoms
            $resn = "SWM4";
            $name =~ s/OW/OH2/;
            $name =~ s/HW1/H1/;
            $name =~ s/HW2/H2/;
            # keep track of how many waters there are
            if ($name =~ /OH2/)
            {
                $nsolv++;
            }
        }
        else
        {
            # it's a SILCS fragment, use resname
            $segid = $resn;
        }

        # check for change in segid here - if it has changed, we reinitialize the residue counter to 1
        if ($atoms > 1)
        {
            my ($prevatom, $prevanum, $prevname, $prevresn, $prevchain, $prevresnr, $prevx, $prevy, $prevz, $prevb, $prevocc) = unpack("A6A5A5A5A1A5A12A8A18A6A6", $in[$i-1]);
            my $prevsegid = '';

            if ( ($prevresn =~ /ADE/) || ($prevresn =~ /CYT/) || ($prevresn =~ /GUA/) || ($prevresn =~ /THY/) || ($prevresn =~ /DA/) || ($prevresn =~ /DC/) || ($prevresn =~ /DG/) || ($prevresn =~ /DT/))
            {
                $prevsegid = "DNAA";
            }
            elsif ($prevresn =~ /POT/)
            {
                $prevsegid = "HETA";
            }
            elsif ($prevresn =~ /SOL/)
            {
                $prevsegid = "SOLV";
            }
            elsif ($prevresn =~ /SWM4/)
            {
                $prevsegid = "SOLV";
            }
            else
            {
                if ($prevresn =~ /IMIA/)
                {
                    $prevresn = "IMID";
                    $prevsegid = "IMID";
                }

                $prevsegid = $prevresn;
                # increment counter for SILCS fragments if a new residue is detected
                if (!($prevsegid =~ $segid))
                {
                    # incr by residue name...
                    if ($resn =~ /ACEY/)
                    {
                        $nacey++;
                    }
                    elsif ($resn =~ /BENX/)
                    {
                        $nbenx++;
                    }
                    elsif ($resn =~ /DMEE/)
                    {
                        $ndmee++;
                    }
                    elsif ($resn =~ /FORM/)
                    {
                        $nform++;
                    }
                    elsif ($resn =~ /IMID/)
                    {
                        $nimid++;
                    }
                    elsif ($resn =~ /MAMY/)
                    {
                        $nmamy++;
                    }
                    elsif ($resn =~ /MEOH/)
                    {
                        $nmeoh++;
                    }
                    elsif ($resn =~ /PRPX/)
                    {
                        $nprpx++;
                    }
                    elsif ($resn =~ /SWM4/)
                    {
                        # do nothing
                    }
                    else
                    {
                        die "Unknown SILCS fragment found: $resn\n";
                    }
                }
                elsif (($prevsegid =~ $segid) && ($prevresnr != $resnr))
                {
                    # incr by residue name...
                    if ($resn =~ /ACEY/)
                    {
                        $nacey++;
                    }
                    elsif ($resn =~ /BENX/)
                    {
                        $nbenx++;
                    }
                    elsif ($resn =~ /DMEE/)
                    {
                        $ndmee++;
                    }
                    elsif ($resn =~ /FORM/)
                    {
                        $nform++;
                    }
                    elsif ($resn =~ /IMID/)
                    {
                        $nimid++;
                    }
                    elsif ($resn =~ /MAMY/)
                    {
                        $nmamy++;
                    }
                    elsif ($resn =~ /MEOH/)
                    {
                        $nmeoh++;
                    }
                    elsif ($resn =~ /PRPX/)
                    {
                        $nprpx++;
                    }
                    elsif ($resn =~ /SWM4/)
                    {
                        # do nothing
                    }
                    else
                    {
                        die "Unknown SILCS fragment found: $resn\n";
                    }
                }
            }

            if (!($prevsegid =~ $segid))
            {
                # we have changed segments, start numbering over
                $offset = $resnr - 1;
            }
        }
        $resnr -= $offset;  # safe because the first segment will always have offset = 0, all subsequent chains have offset != 0

        # element name for printing last field
        my @aname = split('', $name);
        my $elem = $aname[0];
        printf OUT $pdbfmt, $atom, $anum, $name, $resn, $chain, $resnr, $x, $y, $z, $b, $occ, $segid, $elem;
    }
    else
    {
        print OUT $in[$i];
    }
}
print OUT "TER\n";
print OUT "ENDMDL\n";
close(OUT);

# print to stream file
print STROUT "set boxx  $boxx\n";
print STROUT "set boxy  $boxy\n";
print STROUT "set boxz  $boxz\n";
print STROUT "set alp   $alp\n";
print STROUT "set bet   $bet\n";
print STROUT "set gam   $gam\n";
print STROUT "set shape $shape\n\n";
print STROUT "set seq   $seq\n";
print STROUT "set ndnaa $ndnaa\n";
print STROUT "set nheta $nheta\n";
print STROUT "set nsolv $nsolv\n";
print STROUT "set nacey $nacey\n";
print STROUT "set nbenx $nbenx\n";
print STROUT "set ndmee $ndmee\n";
print STROUT "set nform $nform\n";
print STROUT "set nimid $nimid\n";
print STROUT "set nmamy $nmamy\n";
print STROUT "set nmeoh $nmeoh\n";
print STROUT "set nprpx $nprpx\n\n";

close(STROUT);

exit;
