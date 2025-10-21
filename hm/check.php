<?php

chdir(__DIR__);
chdir("..");

require_once("data/protutils.php");
$dofix = @$argv[1] == "fix";
$start = @$argv[2];
$nochk = (@$argv[2] == "nochk" || @$argv[3] == "nochk") ? "nochk" : "";
$errors = 0;
$byfam = [];
$ttlfam = [];

foreach ($prots as $rcpid => $p)
{
    if ($start)
    {
        if (substr($start, 0, 2) == 'OR')
        {
            if (intval(substr($rcpid, 2, 2)) < intval(substr($start, 2, 2))) continue;
            else if (intval(substr($rcpid, 2, 2)) == intval(substr($start, 2, 2)) && $rcpid < $start) continue;
        }
        else if ($rcpid < $start) continue;
    }

    $fam = family_from_protid($rcpid);
    $famno = (substr($fam, 0, 2) == "OR") ? intval(substr($fam, 2)) : 0;
    if ($famno)
    {
        if (!isset($byfam[$famno])) $byfam[$famno] = 0;
        if (!isset($ttlfam[$famno])) $ttlfam[$famno] = 1;
        else $ttlfam[$famno]++;
    }

    if (!file_exists("pdbs/$fam/$rcpid.active.pdb"))
    {
        echo "$rcpid active PDB is missing.\n";
        if ($famno) $byfam[$famno]++;
        $errors++;
        if (file_exists("pdbs/$fam/$rcpid.inactive.pdb") && $dofix) passthru("php -f hm/dohm.php $rcpid $nochk");
        continue;
    }
}

echo "$errors errors.\n";

foreach ($byfam as $famno => $bf)
{
    echo "OR$famno: $bf/{$ttlfam[$famno]} missing.\n";
}
