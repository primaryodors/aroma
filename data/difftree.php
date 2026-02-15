<?php

chdir(__DIR__);
require_once("protutils.php");

$partial_match =
[
    "A" => "GMILV",
    "R" => "KH",
    "N" => "Q",
    "D" => "E",
    "C" => "ST",
    "E" => "Q",
    "Q" => "N",
    "G" => "AP",
    "H" => "RKFWY",
    "I" => "MALV",
    "L" => "MAIV",
    "K" => "RH",
    "M" => "AILV",
    "F" => "HWY",
    "P" => "G",
    "S" => "CT",
    "T" => "ST",
    "W" => "FHY",
    "Y" => "WHF",
    "V" => "MAIL"
];

function compare_aligneds($a, $b)
{
    $alen = strlen($a);
    $blen = strlen($b);
    $len = min($alen, $blen);
    $e = 0;
    for ($i=0; $i<$len; $i++)
    {
        $c = substr($a, $i, 1);
        $d = substr($b, $i, 1);
        if ($c == '-' || $d == '-') continue;
        if ($c != $d) $e++;
    }

    return $e;
}

function make_consensus($inarr)
{
    global $partial_match;
    $result = "";
    $n = strlen($inarr[array_keys($inarr)[0]]);
    for ($i=0; $i<$n; $i++)
    {
        $bins = [];
        $c = "-";
        foreach ($inarr as $str)
        {
            $lc = substr($str, $i, 1);
            if (!isset($bins[$lc])) $bins[$lc] = 1.0;
            else $bins[$lc] += 1.0;

            if (isset($partial_match[$lc]))
            {
                $m = strlen($partial_match[$lc]);
                for ($j=0; $j<$m; $j++)
                {
                    $llc = substr($partial_match[$lc], $j, 1);
                    $pm = 0.25*amino_similarity($lc, $llc);
                    if (!isset($bins[$llc])) $bins[$llc] = $pm;
                    else $bins[$llc] += $pm;
                }
            }
        }
        arsort($bins);
        if (count($bins)) $c = array_keys($bins)[0];
        $result .= $c;
    }

    return $result;
}


$md1 = $md2 = false;
$mdelta = 0;
$alis = [];
foreach ($prots as $rcpid1 => $p1)
{
    if (substr($rcpid1, 0, 2) == "VN") continue;
    if (substr($rcpid1, 0, 2) == "MS") continue;
    $alis[$rcpid1] = $p1['aligned'];
    foreach ($prots as $rcpid2 => $p2)
    {
        if (substr($rcpid2, 0, 2) == "VN") continue;
        if (substr($rcpid2, 0, 2) == "MS") continue;

        // echo "Comparing $rcpid1 ~ $rcpid2...\n";
        $delta = compare_aligneds($p1['aligned'], $p2['aligned']);
        if ($delta > $mdelta)
        {
            $md1 = $rcpid1;
            $md2 = $rcpid2;
            $mdelta = $delta;
            echo "Delta of $md1 ~ $md2: $mdelta\n";
        }
    }
}

echo "Greatest difference: $md1 ~ $md2\n";

$consensus = make_consensus($alis);
echo "Consensus is:\n$consensus\n\n";

$btree = ["b0" => $md1, "b1" => $md2];
$unused = [];
foreach (array_keys($prots) as $rcpid)
{
    if (substr($rcpid, 0, 2) == "VN") continue;
    if (substr($rcpid, 0, 2) == "MS") continue;
    if ($rcpid != $md1 && $rcpid != $md2) $unused[$rcpid] = $rcpid;
}

while (count($unused))
{
    $nearest = false;
    $nnode = false;
    $ndelta = 999999999;
    foreach ($unused as $rcpid)
    {
        $p1 = $prots[$rcpid];
        foreach ($btree as $node => $cmpid)
        {
            if (strlen($cmpid) > 100) $seq2 = $cmpid;
            else
            {
                $p2 = $prots[$cmpid];
                $seq2 = $p2['aligned'];
            }
            $delta = compare_aligneds($p1['aligned'], $seq2);
            // echo "Comparing $rcpid ~ $cmpid delta $delta/$ndelta...\n";
            if ($delta < $ndelta)
            {
                $nearest = $rcpid;
                $nnode = $node;
                $ndelta = $delta;
            }
        }
    }

    if (false===$nearest || false===$nnode) die("Something went wrong.\n");

    $lnn = strlen($nnode);
    $oldn = [];
    $newbt = [];
    foreach ($btree as $lnode => $v)
    {
        if ($lnode != $nnode && substr($lnode, 0, $lnn) == $nnode)
        {
            $oldn[] = $lnode;
            $newbt[$nnode.'0'.substr($lnode, $lnn)] = $v;
        }
    }
    foreach ($oldn as $ln) unset($btree[$ln]);
    foreach ($newbt as $nn => $nv) $btree[$nn] = $nv;

    $btree["{$nnode}0"] = $btree[$nnode];
    $btree["{$nnode}1"] = $nearest;

    // Old node becomes the consensus of the two neighbors
    $base = [$consensus, strlen($btree[$nnode]) > 100 ? $btree[$nnode] : $prots[$btree[$nnode]]['aligned'], $prots[$nearest]['aligned']];
    $btree[$nnode] = make_consensus($base);

    unset($unused[$nearest]);
    echo '.';
    // ksort($btree); print_r($btree);
    // if (count($btree) >= 20) break;
}
echo "\n";

foreach ($btree as $node => $value) if (strlen($value) > 100) unset($btree[$node]);

ksort($btree);
print_r($btree);