<?php

chdir(__DIR__);
require_once("data/protutils.php");
chdir(__DIR__);

foreach ($prots as $protid => $p)
{
    if (substr($protid, 0, 2) == "OR")
    {
        passthru("php -f hm/dohm.php $protid");
    }
}