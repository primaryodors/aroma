<?php

chdir(__DIR__);
require_once("../data/protutils.php");
chdir(__DIR__);

$exp = file_get_contents("experimental.ali"); //."\n\n".file_get_contents("hm.ali");

$deletions =
[
    "5cxv" => "A|1000|1200",
];

$prevln = "";
$already = [];
foreach (explode("\n", $exp) as $ln)
{
    if (substr($prevln, 0, 4) == ">P1;")
    {
        $rcpid = trim(substr($prevln, 4));
        $already[$rcpid] = $rcpid;
        $pettias = explode(':', $ln);
        $pseqid = $pettias[1];
        $fn = "$pseqid.pdb";
        if (!file_exists($fn))
        {
            if (preg_match("/^[0-9A-Za-z]{4}$/", $pseqid))
            {
                $rcsbid = strtoupper($pseqid);
                $url = "https://files.rcsb.org/download/$rcsbid.pdb";
                $c = file_get_contents($url);
                if (isset($deletions[$pseqid]))
                {
                    list($delstrand, $delsr, $deler) = explode('|',$deletions[$pseqid]);
                    $c = explode("\n", $c);
                    foreach ($c as $k => $ln)
                    {
                        if (substr($ln, 0, 6) == "ATOM  ")
                        {
                            $strand = substr($ln, 21, 1);
                            $resno = intval(substr($ln, 22, 4));
                            if ($strand == $delstrand && $resno >= $delsr && $resno <= $deler)
                            {
                                echo "Delete: $ln\n";
                                unset($c[$k]);
                            }
                        }
                    }
                    $c = implode("\n", $c);
                }
                file_put_contents($fn, $c);
                if (file_exists($fn)) echo "Downloaded $rcsbid.\n";
            }
        }
    }

    $prevln = $ln;
}

$fp = fopen("allgpcr.ali", "w");
if (!$fp) die("FAIL; check folder permissions.\n");
fwrite($fp, $exp);
fwrite($fp, "\n\n");

foreach ($prots as $rcpid => $p)
{
    if (!isset($p['aligned'])) continue;
    if (@$already[$rcpid]) continue;

    $paligned = "";
    $temp = $p['aligned'];
    if (substr($temp, -4) == '----') $temp = substr($temp, 0, -4);
    while ($temp)
    {
        $paligned .= substr($temp, 0, 130)."\n";
        $temp = substr($temp, 130);
    }

    $p1row = ">P1;$rcpid";
    fwrite($fp, "$p1row\n");

    $fam = family_from_protid($rcpid);
    $mem = member_from_protid($rcpid);
    $pname = "(unspecified)";
    $gpseq = "";
    switch ($fam)
    {
        case "TAAR":
        $pname = "Trace amine-associated receptor $mem";
        break;

        case "VN1R":
        $pname = "Vomeronasal type 1 receptor number $mem";
        break;

        case "MS4A":
        $pname = "Membrane-spanning 4A receptor $mem";
        break;

        default:
        $famn = substr($fam, 2);
        $sub = subfamily_from_protid($rcpid);
        $pname = "Olfactory receptor family $famn subfamily $sub number $mem";
    }

    $seqlen = strlen($p['sequence']);

    $deets = "sequence:$rcpid:1     :A:$seqlen  :A:$pname:Homo sapiens: 1.90: 0.19";
    fwrite($fp, "$deets\n");
    fwrite($fp, "$paligned---------------------------------------------------------------------------------------------------*\n\n");

    /* $p1row = ">P1;$rcpid.i";
    fwrite($fp, "$p1row\n");
    fwrite($fp, "$deets\n");
    fwrite($fp, "$paligned---------------------------------------------------------------------------------------------------*\n\n"); */
}

fclose($fp);
echo "Wrote alignments file.\n";
