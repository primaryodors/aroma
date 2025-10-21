<?php

chdir(__DIR__);
require_once("odorutils.php");

function num_heavy_atoms($smiles)
{
    return strlen(preg_replace("/[^A-Za-z]/", "", $smiles));
}

function get_canonical($o)
{
    $result = [];
    exec("obabel -:'{$o['smiles']}' -ocan 2>/dev/null", $result);
    if (!count($result)) echo "ERROR no output for {$o['smiles']}\n";
    return trim($result[0]);
}

foreach ($odors as $i => $o1)
{
    $canonical = get_canonical($o1);
    $atoms1 = num_heavy_atoms($canonical);
    foreach ($odors as $j => $o2)
    {
        if ($i == $j) continue;
        $atoms2 = num_heavy_atoms($o2['smiles']);
        if ($atoms2 != $atoms1) continue;

        $canon1 = get_canonical($o2);
        if ($canon1 == $canonical)
        {
            echo $o1['full_name']." and ".$o2['full_name']." may be the same molecule. The canonical SMILES is $canonical.\n";
        }
    }
}