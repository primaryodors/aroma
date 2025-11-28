<?php

chdir(__DIR__);
require_once("../data/protutils.php");
chdir(__DIR__);

if (!file_exists('/usr/lib/libmodeller.so')) die("This feature requires MODELLER.\nPlease see: https://salilab.org/modeller/ to obtain this third party software.\n\n");

$rcpid = @$argv[1];
$p = $prots[$rcpid];
if (!$rcpid) die("Usage:\nphp -f hm/dohm.php PROTID\n\n");

$riglig = $hetatms = false;
if (in_array("riglig", $argv)) $riglig = true;

$disulfs = $hb = $hetatmln = $restraint345 = $helices = "";
foreach ($p["region"] as $rgname => $rgnse)
{
    $nmsub3 = substr($rgname, 0, 3);
    if ($nmsub3 != "TMR" && $nmsub3 != "HXR") continue;
    if (($rcpid == "OR10Q1" || $rcpid == "OR52A1" || $rcpid == "OR52A4") && $nmsub3 == "HXR") continue;
    $rgs = $rgnse['start'];
    $rge = $rgnse['end'];
    $helices .= "        rsr.add(secondary_structure.Alpha(self.residue_range('$rgs:A', '$rge:A')))\n";
}

if (substr($rcpid, 0, 2) == "OR")
{
    $rgs = resno_from_bw($rcpid, "45.52");
    $rge = resno_from_bw($rcpid, "45.58");
    $helices .= "        rsr.add(secondary_structure.Alpha(self.residue_range('$rgs:A', '$rge:A')))\n";
}

$xlinx = [ ["3.25", "45.50"] ];
foreach ($xlinx as $xl)
{
    try
    {
        $rno1 = resno_from_bw($rcpid, $xl[0]);
        $rno2 = resno_from_bw($rcpid, $xl[1]);
    }
    catch (Exception $e)
    {
        continue;
    }

    if (!$rno1 || !$rno2) continue;

    $raa1 = substr($p['sequence'], $rno1-1, 1);
    $raa2 = substr($p['sequence'], $rno2-1, 1);
    if ($raa1 == 'C' && $raa2 == 'C')
    {
        $disulfs .=
            "        self.patch(residue_type='DISU', residues=(self.residues['$rno1:A'],\n"
            ."                                                  self.residues['$rno2:A']))\n";
        $restraint345 .= "        rsr.add(forms.Gaussian(group=physical.xy_distance,
                                feature=features.Distance(at['SG:$rno1:A'],
                                                          at['SG:$rno2:A']),
                                mean=2.05, stdev=0.2))";
    }
}

if ($disulfs) $disulfs = "    def special_patches(self, aln):\n$disulfs";

$mdlcls = "DOPEHRLoopModel";
if ($rcpid == "OR51F1" || $rcpid == "OR51F2" || $rcpid == "OR51G2") $mdlcls = "AutoModel";

exec("php -f build_alignment_file.php");

$consOR51 = "'8uxv'";
$consOR52 = "'8hti'";
$OR51E2 = "'8f76'";
$consOR1 = "'8uxy'";
$consOR2 = "'8uy0'";
$consOR2b = "'8uy0b'";
$consOR4 = "'8uyq'";
$CLASSII = "$consOR1, $consOR2, $consOR4";
$TAAR1 = "'8jln', '8jlo', '8jlp', '8jlq', '8jlr', '8jso'";
$mTAAR = "'8iwe', '8iwm', '8itf', '8iw4', '8iw9', '8pm2'";
$ADORA2A = "'6gdg'";
$ADRB2 = "'7dhr', '8gej'";
$LPAR1 = "'7td0', '7yu3'";
$TAS2R = "'7xp6'";
$CB = "'5xr8', '5xra'";
$CHRM1 = "'6oij'";

$hmrcpstrid = 'A';
$adjustments = "";
$mtl539 = $mtl546 = false;

$mdlevel = "";
if ($rcpid == "OR2AE1" || $rcpid == "OR2AG1" || $rcpid == "OR2AG2") $mdlevel = "a.md_level = refine.very_slow";

$restraints_misc = [];

$fam = family_from_protid($rcpid);
$sub = subfamily_from_protid($rcpid);
$famsub = "$fam$sub";

switch ($fam)
{
    case 'TAAR':
    if ($rcpid == "TAAR1")
    {
        $knowns = $TAAR1;
        $mdlcls = "AutoModel";
    }
    else if ($rcpid == "TAAR9") $knowns = $mTAAR;
    else $knowns = "$mTAAR"; //, $TAAR1";
    break;

    case 'VN1R':
    $knowns = "$TAS2R, $CB";
    break;

    case 'MS4A':
    die("No HM templates known for MS4A receptors.\n");
    break;

    case 'OR1':
    $knowns = "$consOR1";
    break;

    case 'OR2':
    if ($rcpid == "OR2AP1") $knowns = "$consOR4";
    else if ($famsub == "OR2M")
    {
        $knowns = "$CHRM1, $CLASSII";
        $restraints_misc[] = "3.33|5.42|7.3";
    }
    else $knowns = "$CLASSII";
    break;

    case 'OR3':
    $knowns = "$consOR1";
    break;

    case 'OR4':
    $knowns = "$consOR4";
    break;

    case 'OR5':
    if ($rcpid == "OR5V1") $knowns = "$consOR4";
    else $knowns = "$consOR1";
    break;

    case 'OR6':
    $knowns = "$consOR1, $consOR4";
    if ($famsub == "OR6A" || $famsub == "OR6B" || $famsub == "OR6P" || $famsub == "OR6Y")
        $restraints_misc[] = "4.60:NZ|5.39:OD2|2.53";
    break;

    case 'OR7':
    $knowns = "$consOR1";
    break;

    case 'OR8':
    if ($rcpid == "OR8S1") $knowns = "$consOR4";
    else $knowns = "$consOR1";
    break;

    case 'OR9':
    $knowns = "$consOR1";
    break;

    case 'OR10':
    if ($rcpid == "OR10AD1") $knowns = "$consOR2";
    else $knowns = "$consOR4";
    break;

    case 'OR11':
    $knowns = "$consOR4";
    break;

    case 'OR12':
    $knowns = "$consOR4";
    break;

    case 'OR13':
    if ($famsub == "OR13A" || $famsub == "OR13G") $knowns = "$consOR1";             // OR13A/G are actually OR3s.
    else $knowns = "$consOR2";
    break;

    case 'OR14':
    $knowns = "$consOR2";
    break;

    case 'OR51':
    if ($rcpid == "OR51E2")
    {
        $knowns = false;
        $pyoutfn = "8f76.pdb";
        $hmrcpstrid = 'A';
    }
    else $knowns = "$consOR51";
    break;

    case 'OR52':
    $knowns = "$consOR52";
    break;

    case 'OR56':
    $knowns = "$consOR51, $consOR52";
    break;

    default:
    $knowns = "$consOR1, $consOR2, $consOR4, $consOR51, $consOR52, $mTAAR, $TAAR1";
}

if ($knowns)
{
    $atom655 = $atom4551 = false;
    $restraint456 = "";
    if (substr($rcpid, 0, 2) == "OR" && intval(substr($fam, 2)) < 50)
    {
        $rno655 = resno_from_bw($rcpid, "6.55");
        $aa655 = substr($p['sequence'], $rno655-1, 1);
        if ($aa655 == 'Y') $atom655 = "OH:$rno655:A";
        else if ($aa655 == 'H') $atom655 = "NE2:$rno655:A";

        $rno4551 = resno_from_bw($rcpid, "45.51");
        $rno4552 = $rno4551 + 1;
        $aa4551 = substr($p['sequence'], $rno4551-1, 1);
        $aa4552 = substr($p['sequence'], $rno4552-1, 1);
        if ($aa4551 == 'D') $atom4551 = "OD1:$rno4551:A";
        else if ($aa4551 == 'N') $atom4551 = "OD1:$rno4551:A";
        else if ($aa4551 == 'E') $atom4551 = "OE1:$rno4551:A";
        else if ($aa4551 == 'Q') $atom4551 = "OE1:$rno4551:A";
        else if ($aa4552 == 'D') $atom4551 = "OD1:$rno4552:A";
        else if ($aa4552 == 'E') $atom4551 = "OE1:$rno4552:A";

        if ($aa4551 == 'E' || $aa4551 == 'Q' || $aa4552 == 'E') $knowns = "$consOR2b";
    }

    if ($atom655 && $atom4551)
        $restraint456 = "        rsr.add(forms.Gaussian(group=physical.xy_distance,
                                feature=features.Distance(at['$atom655'],
                                                          at['$atom4551']),
                                mean=2.5, stdev=0.5))";

    $atom558 = $atom753 = false;
    $restraint57 = "";

    $rno558 = resno_from_bw($rcpid, "5.58");
    $aa558 = substr($p['sequence'], $rno558-1, 1);
    if ($aa558 == 'Y') $atom558 = "OH:$rno558:A";

    $rno753 = resno_from_bw($rcpid, "7.53");
    $aa753 = substr($p['sequence'], $rno753-1, 1);
    if ($aa753 == 'Y') $atom753 = "OH:$rno753:A";

    if ($atom558 && $atom753)
        $restraint57 = "        rsr.add(forms.Gaussian(group=physical.xy_distance,
                                feature=features.Distance(at['$atom558'],
                                                          at['$atom753']),
                                mean=4.6, stdev=1.2))";

    $restraintmtl = "";
    if ($famsub == "OR2M" || $famsub == "OR2T")
    {
        $rno539 = resno_from_bw($rcpid, "5.39");
        $rno542 = resno_from_bw($rcpid, "5.42");
        $rno543 = resno_from_bw($rcpid, "5.43");
        $rno546 = resno_from_bw($rcpid, "5.46");

        $aa539 = substr($p['sequence'], $rno539-1, 1);
        $aa542 = substr($p['sequence'], $rno542-1, 1);
        $aa543 = substr($p['sequence'], $rno543-1, 1);
        $aa546 = substr($p['sequence'], $rno546-1, 1);

        if ($aa542 == 'C' && $aa543 == 'C')
        {
            $cmdist = 4.7;
            $atom542 = "SG:$rno542:A";
            $atom543 = "SG:$rno543:A";
            if ($aa539 == 'M')
            {
                $mtl539 = true;
                $atom539 = "SD:$rno539:A";
                $restraintmtl .= "        rsr.add(forms.Gaussian(group=physical.xy_distance,
                                feature=features.Distance(at['$atom539'],
                                                            at['$atom542']),
                                mean=$cmdist, stdev=0.25))\n\n";
                $restraintmtl .= "        rsr.add(forms.Gaussian(group=physical.xy_distance,
                                feature=features.Distance(at['$atom539'],
                                                            at['$atom543']),
                                mean=$cmdist, stdev=0.25))\n\n";
            }
            if ($aa546 == 'M')
            {
                $mtl546 = true;
                $atom546 = "SD:$rno546:A";
                $restraintmtl .= "        rsr.add(forms.Gaussian(group=physical.xy_distance,
                                feature=features.Distance(at['$atom546'],
                                                            at['$atom542']),
                                mean=$cmdist, stdev=0.25))\n\n";
                $restraintmtl .= "        rsr.add(forms.Gaussian(group=physical.xy_distance,
                                feature=features.Distance(at['$atom546'],
                                                            at['$atom543']),
                                mean=$cmdist, stdev=0.25))\n\n";
            }
            if ($mtl539 || $mtl546)
            {
                $restraintmtl .= "        rsr.add(forms.Gaussian(group=physical.xy_distance,
                                feature=features.Distance(at['$atom542'],
                                                            at['$atom543']),
                                mean=3.8, stdev=0.2))\n\n";
            }
        }
    }

    $restraints_misc_str = "";
    foreach ($restraints_misc as $rm)
    {
        list($bw1, $bw2, $rmr) = explode('|',$rm);
        if (false!==strpos($bw1, ":")) list($bw1, $a1) = explode(':', $bw1);
        else $a1 = "CA";
        if (false!==strpos($bw2, ":")) list($bw2, $a2) = explode(':', $bw2);
        else $a2 = "CA";
        $rno1 = resno_from_bw($rcpid, $bw1);
        $rno2 = resno_from_bw($rcpid, $bw2);
        if (!$rno1 || !$rno2) continue;
        $rmr = floatval($rmr);
        $restraints_misc_str .= "        rsr.add(forms.Gaussian(group=physical.xy_distance,
                                feature=features.Distance(at['$a1:$rno1:A'],
                                                          at['$a2:$rno2:A']),
                                mean=$rmr, stdev=1.2))\n";
    }

    if ($riglig)
    {
        if (file_exists("riglig.pdb")) unlink("riglig.pdb");
        $rigpdb = "../pdbs/$fam/$rcpid.riglig.pdb";
        if (!file_exists($rigpdb)) die("ERROR: No riglig file for $rcpid.\n\nTo create a riglig file:\n1.) Run the prediction for a known agonist in a receptor, allowing the prediction to fail, or you can interrupt it when the splash graphic appears;\n2.) You will see the exact AromaDock command-line call used, located immediately above the splash graphic;\n3.) Copy-paste this command (it will begin with \"bin/aromadock\") and add the --movie argument, then run it and allow it to finish;\n4.) Using the 3D viewer (viewer.htm) in a web browser, click \"Choose File\" at upper left;\n5.) Navigate to the aromadock/tmp folder and open _iters.dock;\n6.) Look through all poses and nodes to find the best ligand positioning, either by visual examination or by the energies given by clicking the Stats (".mb_chr(0x1f4ca).") icon;\n7.) Remembering the pose and node numbers of the best positioning, open the tmp/_iters.dock file in a text editor and find the right section of the file.\n\tFor example, if the desired pose is number 10 and the node is number 3, look for Pose: 10 followed by Node: 3 in the _iters.dock.\n\tHint: It will probably be easier to skip ahead to the next section, e.g. Pose: 10 Node: 4 in this example, and scroll up.\n8.) Find the collection of lines that start with HETATM; this is your best-positioned ligand;\n9.) Copy-paste the HETATM lines to a new text file and replace all occurrences of \"LIG     0\" (five spaces between LIG and 0) with \"LIG A   1\" (three spaces between A and 1);\n10.) Add a newline followed by TER then another newline then END and one last newline at the end of the text file and save it in aromadock/pdbs/{family}/{receptor}.riglig.pdb\n\tFor example, if this is for OR1A1 then you would save the riglig file as aromadock/pdbs/OR1/OR1A1.riglig.pdb\n11.) Rerun this script for that receptor with the riglig argument.\n\n");
        copy($rigpdb, "riglig.pdb");
        $knowns .= ", 'riglig'";
        $hetatms = true;
    }

    if ($hetatms) $hetatmln = "env.io.hetatm = True";

    // Generate Python script
    $py = <<<natrixs

from modeller import *
from modeller.automodel import *

env = Environ()

$hetatmln

class MyModel($mdlcls):
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms
$helices
$restraint345
$restraint456
$restraint57
$restraintmtl
$restraints_misc_str
$disulfs

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

a = MyModel(env,
              alnfile  = 'allgpcr.ali',
              knowns   = ($knowns),
              sequence = '$rcpid')
a.starting_model = 0
a.ending_model   = 9
a.library_schedule = autosched.slow
a.max_var_iterations = 300
$mdlevel
a.make()

natrixs;

    if (file_exists("$rcpid.hm.py"))
    {
        $result = [];
        exec("ps -ef | grep -E \"python[0-9] $rcpid.hm.py\" | grep -v grep", $result);
        if (count($result)) die("Already running.\n");
    }

    $fp = fopen("$rcpid.hm.py", "w");
    fwrite($fp, $py);
    fclose($fp);

    @unlink("hm.out");
    passthru("python3 $rcpid.hm.py | tee $rcpid.hm.out");
    $c = file_get_contents("$rcpid.hm.out");
    $best_energy = 1e9;
    $pyoutfn = false;
    $mode = false;
    foreach (explode("\n", $c) as $ln)
    {
        if (false !== strpos($ln, "Summary of successfully produced models:")) $mode = true;
        if (false !== strpos($ln, "Summary of successfully produced loop models:")) $mode = true;
        if (preg_match("/Filename\\s+molpdf/i", $ln)) $mode = true;
        if (!$mode) continue;

        $ln = preg_replace("/\\s+/", " ", $ln);
        $pieces = explode(" ", $ln);
        if (count($pieces) < 2) continue;
        if (!is_numeric($pieces[1])) continue;
        $e = floatval($pieces[1]);
        if ($e < $best_energy)
        {
            $best_energy = $e;
            $pyoutfn = $pieces[0];
        }
    }

    if (!$pyoutfn) die("FAIL.\n");

    $famno = intval(preg_replace("/[^0-9]/", "", $fam));
    $adjustments = "";
    if ($famno < 50) $adjustments .= "IF $3.37 != \"G\" THEN ATOMTO %3.37 EXTENT @6.48\n";
    else if ($famno == 51 || $famno == 52) $adjustments .= "ATOMTO %6.59 EXTENT @4.57\n";
    else if ($rcpid == "OR56B2") $adjustments .= "ATOMTO %6.58 EXTENT @4.57\n";

    $knowns = preg_replace("/[^0-9a-zA-Z_ ]/", "", $knowns);
}

$phew = <<<blixtos

LET \$rcpid = "$rcpid"

LET \$inpf = "pdbs/$fam/$rcpid.inactive.pdb"

LET \$mdld = "hm/$pyoutfn"

LOAD \$inpf A I
LET %rcpseqln = %SEQLENI
LOAD \$mdld $hmrcpstrid A

BWCOPY I A

STRAND I
UPRIGHT I
BWCENTER

STRAND A

REMARK   1
REMARK   1 REFERENCE 1
REMARK   1  AUTH 1 B. Webb, A. Sali.
REMARK   1  TITL 1 Comparative Protein Structure Modeling Using Modeller.
REMARK   1  REF  1 Current Protocols in Bioinformatics 54, John Wiley & Sons, Inc., 5.6.1-5.6.37, 2016.
REMARK   1  DOI  1 10.1002/0471250953.bi0506s15
REMARK   1  
REMARK   1  AUTH 2 M.A. Marti-Renom, A. Stuart, A. Fiser, R. Sánchez, F. Melo, A. Sali.
REMARK   1  TITL 2 Comparative protein structure modeling of genes and genomes.
REMARK   1  REF  2 Annu. Rev. Biophys. Biomol. Struct. 29, 291-325, 2000.
REMARK   1  DOI  2 10.1146/annurev.biophys.29.1.291
REMARK   1
REMARK   1  AUTH 3 A. Sali & T.L. Blundell.
REMARK   1  TITL 3 Comparative protein modelling by satisfaction of spatial restraints.
REMARK   1  REF  3 J. Mol. Biol. 234, 779-815, 1993.
REMARK   1  DOI  3 10.1006/jmbi.1993.1626
REMARK   1  
REMARK   1  AUTH 4 A. Fiser, R.K. Do, & A. Sali.
REMARK   1  TITL 4 Modeling of loops in protein structures.
REMARK   1  REF  4 Protein Science 9. 1753-1773, 2000.
REMARK   1  DOI  4 10.1110/ps.9.9.1753
REMARK   1  
IF "$knowns" = "" REMARK 265 HM_TEMPLATES: none
ELSE REMARK 265 HM_TEMPLATES: $knowns

DELETE 1 %1.26
HYDRO

UNCHAIN I
UNCHAIN O
STRAND A
UPRIGHT
BWCENTER
$adjustments

IF $3.25 != "C" OR $45.50 != "C" GOTO _not_disulfide
DELATOM %3.25 HG
DELATOM %45.50 HG
_not_disulfide:

IF $5.42 != "C" OR $5.43 != "C" GOTO _not_Cu_binding_site
IF $5.39 != "M" GOTO _not_Cu_539
MEASURE %5.39 "SD" %5.42 "SG" &sdist
ECHO "5.39:SD - 5.42:SG distance: " &sdist
MEASURE %5.39 "SD" %5.43 "SG" &sdist
ECHO "5.39:SD - 5.43:SG distance: " &sdist
_not_Cu_539:
IF $5.46 != "M" GOTO _not_Cu_539
MEASURE %5.46 "SD" %5.42 "SG" &sdist
ECHO "5.46:SD - 5.42:SG distance: " &sdist
MEASURE %5.46 "SD" %5.43 "SG" &sdist
ECHO "5.46:SD - 5.43:SG distance: " &sdist
_not_Cu_546:
MEASURE %5.42 "SG" %5.43 "SG" &sdist
ECHO "5.42:SD - 5.43:SG distance: " &sdist
_not_Cu_binding_site:

LET \$outf = "pdbs/$fam/$rcpid.active.pdb"
SAVE \$outf

blixtos;

$fp = fopen("$rcpid.hm.phew", "w") or die("FAILED to open $rcpid.hm.phew for writing.");
fwrite($fp, $phew);
fclose($fp);

chdir(__DIR__);
chdir("..");
passthru("./bin/phew hm/$rcpid.hm.phew");

if ($famsub == "OR5K")
{
    passthru("./bin/phew hm/fivewinder.phew pdbs/$fam/$rcpid.active.pdb 4 7");
    passthru("./bin/phew hm/fivewinder.phew pdbs/$fam/$rcpid.active.pdb -3 6");
}

chdir(__DIR__);
// Only delete temporary files if the HM was successful. Few things are more frustrating than finding an unsuccessful HM
// and having to wait for MODELLER to do its thing *again* because the infernal PHP script deleted all the working files.
if (file_exists("../pdbs/$fam/$rcpid.active.pdb") && filemtime("../pdbs/$fam/$rcpid.active.pdb") > filemtime("$rcpid.hm.phew"))
{
    if (!in_array("nodel", $argv)) foreach (glob("$rcpid.*") as $doomed) unlink($doomed);
}
