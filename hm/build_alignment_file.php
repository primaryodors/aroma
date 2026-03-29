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
        $fn = "tpl/$pseqid.pdb";
        if (!file_exists($fn))
        {
            if (preg_match("/^[0-9A-Za-z]{4}$/", $pseqid))
            {
                $rcsbid = strtoupper($pseqid);
                $url = "https://files.rcsb.org/download/$rcsbid.pdb";
                die($url);
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

$cpl = "";
chdir(__DIR__);
if (file_exists("../coupled"))
{
    if (file_exists("../coupled/coupled.ali"))
    {
        $cpl = file_get_contents("../coupled/coupled.ali");
    }
    else
    {
        $lfam = $lsub = "zero";
        foreach ($prots as $rcpid => $p)
        {
            $fam = family_from_protid($rcpid);
            $sub = subfamily_from_protid($rcpid);
            if ($fam == $lfam && $sub == $lsub) continue;

            chdir(__DIR__);
            $path = "../coupled/$fam/$sub";
            if (!file_exists($path)) continue;
            $d = dir($path);
            $files = [];
            while (false !== ($entry = $d->read()))
            {
                if (substr($entry, -4) != ".pdb") continue;
                $pieces = explode('~', substr($entry, 0, -4));
                if (count($pieces) < 2) continue;
                $files[] = $entry;
            }
            natsort($files);
            foreach ($files as $entry)
            {
                echo "Processing $entry...\n";
                $pieces = explode('~', substr($entry, 0, -4));
                list($rcpid, $gprot) = $pieces;

                chdir("..");
                $phew = <<<phew
LOAD "coupled/$fam/$sub/$entry" A A
ECHO \$SEQUENCEA
phew;
                file_put_contents("tmp/sequence.phew", $phew);
                $result = [];
                exec("bin/phew tmp/sequence.phew", $result);
                unlink("tmp/sequence.phew");
                chdir(__DIR__);
                $cplseq = preg_replace("/\\s+/", "", implode("\n", $result));

                // echo "$cplseq\n";

                if (isset($prots[$rcpid]["aligned"]))
                {
                    $append = "";
                    $cpllen = strlen($cplseq);
                    $rcpseq = $prots[$rcpid]["sequence"];
                    $seqlen = strlen($rcpseq);
                    if ($cpllen != $seqlen)
                    {
                        echo "ERROR: SEQUENCE LENGTH MISMATCH $entry:\n$rcpseq\n$cplseq\n\n";
                        continue;
                    }

                    $cplaln = "";
                    $rcpaln = $prots[$rcpid]["aligned"];
                    $n = strlen($rcpaln);
                    $j = 0;
                    for ($i=0; $i<$n; $i++)
                    {
                        $c = substr($rcpaln, $i, 1);
                        if ($c >= 'A' && $c <= 'Z')
                        {
                            $cplaln .= substr($cplseq, $j, 1);
                            $j++;
                        }
                        else $cplaln .= $c;
                    }

                    $caligned = "";
                    $temp = $cplaln;
                    if (substr($temp, -4) == '----') $temp = substr($temp, 0, -4);
                    while ($temp)
                    {
                        $caligned .= substr($temp, 0, 130)."\n";
                        $temp = substr($temp, 130);
                    }
                    $caligned .= $append;

                    $p1row = ">P1;$rcpid~$gprot";
                    $famno = intval(preg_replace("/[^0-9]/", "", $fam));
                    $memno = substr($rcpid, strlen($fam)+strlen($sub));
                    $structrow = "structure:$rcpid~$gprot:FIRST:A:LAST :A:Olfactory receptor family $famno subfamily $sub member $memno:Homo sapiens: 2.00: 0.20";

                    $new = "$p1row\n$structrow\n$caligned---------------------------------------------------------------------------------------------------*\n\n";
                    // die($new);
                    $cpl .= $new;
                }
                else
                {
                    continue;
                    // TODO:
                }
            }

            $lfam = $fam;
            $lsub = $sub;
            // break;
        }
    }
}

chdir(__DIR__);
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

fwrite($fp, "\n\n");
fwrite($fp, $cpl);

fclose($fp);
echo "Wrote alignments file.\n";
