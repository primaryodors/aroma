<?php

chdir(__DIR__);
include_once("protutils.php");

$partial_match =
[
    "A" => "GMILVSTC",
    "R" => "KH",
    "N" => "Q",
    "D" => "E",
    "C" => "AST",
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
    "S" => "CAT",
    "T" => "SAT",
    "W" => "FHY",
    "Y" => "WHF",
    "V" => "MAIL"
];

$bsrI = ["3.33", "4.57", "4.60", "45.52", "45.53", "5.39", "5.43", "6.55", "6.59"];
$bsrII = ["2.53", "3.41", "3.40", "3.37", "3.36", "3.33", "3.32", "3.29", "4.53", "4.57", "4.60", "45.49", "45.51", "45.52", "5.47", "5.46", "5.43", "5.42", "5.39", "6.48", "6.51", "6.55", "7.42", "7.39", "7.38"];
$bsrT = ["3.29", "3.32", "3.33", "3.36", "3.37", "4.56", "4.61", "45.64", "5.43", "6.48", "6.51", "7.39", "7.42", "7.43"];

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
                    $pm = amino_similarity($lc, $llc);
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

$most = 90;
$more = 80;

foreach ($treenodes as $node => $name)
{
    $node = preg_replace("/[^01]/", "", $node);

    if (substr($node, 0, 3) === "001") $bsr = $bsrT;
    else if (substr($node, 0, 4) === "0001") $bsr = $bsrI;
    else if (substr($node, 0, 4) === "0000") $bsr = $bsrII;
    else continue;
    sort($bsr, SORT_STRING);

    $protbsr = [];
    foreach ($prots as $rcpid => $prot)
    {
        if (preg_match("/^$node/", @$prot['btree']))
        {
            $protbsr[$rcpid] = [];
            foreach ($bsr as $bw)
            {
                $aa = letter_at_bw($rcpid, $bw);
                $protbsr[$rcpid][$bw] = $aa;
            }
        }
    }
    $nprots = count($protbsr);
    if ($nprots < 2) continue;

    $bins = [];
    foreach ($bsr as $bw)
    {
        $bins[$bw] = [];
        foreach ($protbsr as $rcpid => $letters)
        {
            $aa = $letters[$bw];
            if (!isset($bins[$bw][$aa])) $bins[$bw][$aa] = 1.0;
            else $bins[$bw][$aa] += 1.0;

            // foreach ($bins[$bw] as $aa2 => $v) $bins[$bw][$aa2] += amino_similarity($aa, $aa2);
        }

        arsort($bins[$bw]);
        // print_r($bins[$bw]);
    }

    $top = "";
    $almost = "";
    $lhxm = $lhxa = 0;
    foreach ($bsr as $bw)
    {
        $aa = array_keys($bins[$bw])[0];
        $pcnt = floor($bins[$bw][$aa] * 100/$nprots);
        $helix = intval(explode('.', $bw)[0]);
        if ($pcnt >= $most)
        {
            $space = ($lhxm && $lhxm!=$helix) ? "\n\t" : " ";
            $top .= ($top?",$space":"") . "$aa$bw ($pcnt%)";
            $lhxm = $helix;
        }
        else if ($pcnt >= $more)
        {
            $space = ($lhxa && $lhxa!=$helix) ? "\n\t" : " ";
            $almost .= ($almost?",$space":"") . "$aa$bw ($pcnt%)";
            $lhxa = $helix;
        }
    }
    echo "Group $name ($nprots proteins):\n";
    if ($most) echo "At least $most%:\n\t$top\n";
    if ($almost) echo "At least $more%:\n\t$almost\n";
    echo "\n";
}
