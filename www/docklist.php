<?php
chdir(__DIR__);
require_once("../data/protutils.php");
require_once("../data/odorutils.php");
require_once("dlmenu.php");

$extra_js = ['js/tabs.js'];
$extra_css = ['assets/tabs.css'];
include("header.php");

chdir(__DIR__);


$fh = 44;
$fw = 30;
$nih = 25;
$nt = 20;
$noh = $nih + $nt;
$nw = 18;
$filter_svgdat = "m 0,0 $fw,$fh 0,$nih $nw,$nt 0,-$noh $fw,-$fh Z"

?>
<style>
<?php output_dlmenu_style(); ?>

.liglist tr.trbk_active
{
    background-color: #341!important;
}

.liglist tr.trbk_inactive
{
    background-color: #235!important;
}
</style>
<script>
<?php output_dlmenu_script(); ?>
</script>
<h1>Completed Docks</h1>

<?php 
// echo "<svg height=\"81px\" viewBox=\"0 0 80 90\" xmlns=\"http://www.w3.org/2000/svg\"><path fill=\"magenta\" d=\"$filter_svgdat\"></path></svg>";
?>

<div class="box">
<div class="row content scrollh">

<?php if (isset($_REQUEST['r']) || isset($_REQUEST['o'])) { ?>
<a href="docklist.php">Clear filters</a>
<?php } ?>

<table class="liglist">
    <tr><th>Receptor</th>
        <th>Odorant</th>
        <th>Mode</th>
        <th>Dock Energies</th>
        <th>Poses</th>
        <th>Agonist?</td>
    </tr>
<?php

$frcp = false;
$flig = false;
chdir(__DIR__);
foreach ($prots as $protid => $p)
{
    $fam = family_from_protid($protid);
    $dockpath = "../output/$fam/$protid";
    if (!file_exists($dockpath)) continue;
    $dir = dir($dockpath);
    $files = [];
    while (false!==($fname=$dir->read())) $files[] = $fname;
    natsort($files);

    foreach ($files as $fname)
    {
        if (substr($fname, -5) != ".dock") continue;
        if (false===strpos($fname, "~")) continue;
        list($odor, $mode, $opfisehciet) = explode('.', explode('~', $fname)[1], 3);
        echo "<tr class=\"trbk_$mode\">\n";

        echo "<td><a href=\"receptor.php?r=$protid\">$protid</a>";
        if ($frcp != $protid)
            echo " <a href=\"docklist.php?r=$protid\"><svg height=\"13px\" viewBox=\"0 0 80 90\" xmlns=\"http://www.w3.org/2000/svg\"><path fill=\"#50cea8\" d=\"$filter_svgdat\"></path></svg></a>";
        echo "</td>\n";
        $frcp = $protid;

        $o = find_odorant($odor);
        $fn = $o['full_name'];
        $fnu = str_replace(' ', '_', $fn);
        echo "<td><a href=\"odorant.php?o={$o['oid']}\">$fn</a>";
        if ($flig != $o['oid'])
            echo " <a href=\"docklist.php?o={$o['oid']}\"><svg height=\"13px\" viewBox=\"0 0 80 90\" xmlns=\"http://www.w3.org/2000/svg\"><path fill=\"#50cea8\" d=\"$filter_svgdat\"></path></svg></a>";
        echo "</td>\n";
        $flig = $o['oid'];

        echo "<td>$mode</td>";

        $c = file_get_contents("../output/$fam/$protid/$fname");
        $lines = explode("\n", $c);
        $benerg = 0;
        foreach ($lines as $ln) if (substr($ln, 0, 7) == "Total: ")
        {
            $benerg = floatval(substr($ln, 7));
            break;
        }
        $nump = 0;
        foreach ($lines as $ln)
        {
            if (substr($ln, 0, 6) == "Pose: ") $nump++;
        }
        if (!$nump) $nump = "-";

        echo "<td><a href=\"viewer.php?view=dock&prot=$protid&odor=$fnu&mode=$mode\" target=\"_dock\">";
        echo round($benerg, 4);
        echo "</td>";

        echo @"<td>$nump</td>\n";

        $agonist = "?";
        $pair = best_empirical_pair($protid, $odor, true);
        if ($pair)
        {
            if (isset($pair["adjusted_curve_top"]))
            {
                $act = floatval($pair["adjusted_curve_top"]);
                if ($act > 0) $agonist = "Y";
                else if ($act < 0) $agonist = "N";
                else $agonist = "inv";
            }
            elseif (isset($pair["ec50"]))
            {
                $ec50 = floatval($pair["ec50"]);
                if ($ec50 < 0) $agonist = "Y";
                else $agonist = "N";
            }
            elseif (isset($pair["type"]))
            {
                switch ($pair["type"])
                {
                    case "vsa": case "sa": case "ma": case "wa": case "vwa":
                        $agonist = "Y";
                        break;

                    case "ia":
                        $agonist = "inv";
                        break;

                    default:
                        $agonist = "N";
                }
            }
            elseif (@$pair["antagonist"] == "Y") $agonist = "ant";
        }
        echo @"<td>$agonist</td>\n";
    }
}

?></table>
<?php

output_dlmenu_div();
