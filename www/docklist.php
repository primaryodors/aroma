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

<?php
$args = [];
foreach ($_REQUEST as $k => $v)
{
    if ($k != 'e') $args[] = "$k=$v";
}
$args = implode("&", $args);
?>

<a href="docklist.php?<?php echo $args; ?>">All</a>
<a href="docklist.php?e=1&<?php echo $args; ?>">Only empirical</a>

<?php if (isset($_REQUEST['r']) || isset($_REQUEST['o'])) { ?>
<a href="docklist.php">Clear filters</a>
<?php } ?>

<table class="liglist">
    <tr><th>Receptor</th>
        <th>Odorant</th>
        <th>Dock Energies</th>
        <th>Occlusion</th>
        <th>Poses</th>
        <th>Agonist?</td>
        <th>Predicted</td>
    </tr>
<?php

$cached = [];
$cachemt = 0;
$cachefn = "docklist.cache.json";
if (file_exists($cachefn))
{
    $cachemt = filemtime($cachefn);
    $cached = json_decode(file_get_contents($cachefn), true);
}

$frcp = false;
$flig = false;
chdir(__DIR__);

$graphdat =
[
    0 => [],
    1 => [],
    2 => [],
    3 => []
];
foreach ($prots as $protid => $p)
{
    if (isset($_REQUEST['r']) && $protid != $_REQUEST['r']) continue;
    $fam = family_from_protid($protid);
    $dockpath = "../output/$fam/$protid";
    if (!file_exists($dockpath)) continue;
    $dir = dir($dockpath);
    $files = [];
    while (false!==($fname=$dir->read())) $files[] = $fname;
    natsort($files);

    $rows = [];
    foreach ($files as $fname)
    {
        if (substr($fname, -5) != ".dock") continue;
        if (false===strpos($fname, "~")) continue;
        list($odor, $mode, $opfisehciet) = explode('.', explode('~', $fname)[1], 3);
        $o = find_odorant($odor);

        if (isset($_REQUEST['o']) && $o['oid'] != $_REQUEST['o']) continue;

        if (@$_REQUEST['e'])
        {
            if (!isset($o['activity'])) continue;
            $found = false;
            foreach ($o['activity'] as $url => $ldat)
            {
                if (isset($ldat[$protid])) $found = true;
                if ($found) break;
            }
            if (!$found) continue;
        }

        $rowid = "$protid~$odor";
        if (!isset($rows[$rowid])) $rows[$rowid] = [];

        $fpn = "../output/$fam/$protid/$fname";

        if (isset($cached[$fname]) && filemtime($fpn) < $cachemt)
        {
            extract($cached[$fname]);
        }
        else
        {
            $c = file_get_contents($fpn);
            $lines = explode("\n", $c);
            $benerg = 0;
            $lsfe = $lsbe = $phfe = $phbe = 0;
            $nump = 0;
            $occl = 0;
            foreach ($lines as $ln) 
            {
                if (!$benerg && substr($ln, 0, 7) == "Total: ")
                {
                    $benerg = floatval(substr($ln, 7));
                }
                else if (!$lsfe && substr($ln, 0, 25) == "Ligand solvation energy: ")
                {
                    $lsfe = floatval(substr($ln, 25));
                }
                else if (!$lsbe && substr($ln, 0, 32) == "Ligand pocket solvation energy: ")
                {
                    $lsbe = floatval(substr($ln, 32));
                }
                else if (!$phfe && substr($ln, 0, 25) == "Pocket hydration energy: ")
                {
                    $phfe = floatval(substr($ln, 25));
                }
                else if (!$phbe && substr($ln, 0, 31) == "Pocket bound hydration energy: ")
                {
                    $phbe = floatval(substr($ln, 31));
                }
                else if (!$occl && substr($ln, 0, 25) == "Ligand pocket occlusion: ")
                {
                    $occl = floatval(substr($ln, 25));
                }
                else if (substr($ln, 0, 6) == "Pose: ") $nump++;
            }
            if (!$nump) $nump = "-";

            $cached[$fname] =
            [
                "benerg" => $benerg,
                "lsfe" => $lsfe,
                "lsbe" => $lsbe,
                "phfe" => $phfe,
                "phbe" => $phbe,
                "occl" => $occl,
                "nump" => $nump,
            ];
        }

        $rows[$rowid]["benerg_$mode"] = $benerg + ($lsbe - $lsfe) + ($phbe - $phfe);
        $rows[$rowid]["nump_$mode"] = $nump;
        $rows[$rowid]["occl_$mode"] = $occl;

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
                    case "a": case "vsa": case "sa": case "ma": case "wa": case "vwa":
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
        $rows[$rowid]["agonist"] = $agonist;
    }

    file_put_contents($cachefn, json_encode_pretty($cached));

    foreach ($rows as $k => $r)
    {
        list($protid, $odor) = explode("~", $k);
        extract($r);
        echo "<tr>\n";

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

        echo "<td>";
        echo "<a href=\"viewer.php?view=dock&prot=$protid&odor=$fnu&mode=active\" target=\"_dock\">";
        echo round($benerg_active, 4) ?: "-";
        echo "</a>";
        echo " / ";
        echo "<a href=\"viewer.php?view=dock&prot=$protid&odor=$fnu&mode=inactive\" target=\"_dock\">";
        echo round($benerg_inactive, 4) ?: "-";
        echo "</a>";
        echo "</td>";

        echo @"<td>" . (round($occl_active, 3) ?: "-") . " / " . (round($occl_inactive, 3) ?: "-") . "</td>\n";
        echo @"<td>" . ($nump_active ?: "-") . " / " . ($nump_inactive ?: "-") . "</td>\n";

        echo @"<td>$agonist</td>\n";

        /*
            let o = occlusion, H = active enthalpy

            if H >= 0:
                    o >= 1 - (15 - H)/100
            if -15 < H < 0:
                    o >= 0.85
            if H <= -15:
                    o >= 0.85 - (H + 15)/100
        */

        $prediction = 0;
        $benerg_active = floatval($benerg_active);
        $occl = floatval($occl);
        if ($benerg_active > 0)             $prediction = max(0, 100.0 * -($occl - (17 - $benerg_active)/100));
        else if ($benerg_active > -15)      $prediction = max(0, 100.0 *  ($occl - 0.85));
        else                                $prediction = max(0, 100.0 * -($occl - ($benerg_active + 17)/100));

        echo "<td>$prediction</td>\n";

        if ($agonist == 'Y')
        {
            $graphdat[0][] = [$benerg_active, $occl_active, $prediction];
            $graphdat[2][] = [$benerg_active - $benerg_inactive, $occl_active, $prediction];
        }
        else if ($agonist == 'N')
        {
            $graphdat[1][] = [$benerg_active, $occl_active, $prediction];
            $graphdat[3][] = [$benerg_active - $benerg_inactive, $occl_active, $prediction];
        }
    }
}

?></table>

<table>
    <tr>
        <th>Active enthalpy vs. occlusion, agonists</th>
        <th>Active enthalpy vs. occlusion, non-agonists</th>
        <th>Active &Delta;H vs. occlusion, agonists</th>
        <th>Active &Delta;H vs. occlusion, non-agonists</th>
    </tr>
    <tr>
        <?php
        foreach ($graphdat as $gi => $ldat)
        {
            echo "<td>\n";
            $grid = [];
            for ($y=0; $y<20; $y++)
            {
                $grid[$y] = [];
                for ($x=0; $x<20; $x++)
                    $grid[$y][$x] = 0;
            }
            $gridmax = 1;

            foreach ($ldat as $d)
            {
                $x = max(0, min(19, intval(10 + $d[0]/5)));
                $y = max(0, min(19, intval(20 - $d[1]*20)));

                $grid[$y][$x] += 1;
                if ($grid[$y][$x] > $gridmax) $gridmax = $grid[$y][$x];
            }

            echo "Peak = $gridmax<br>\n";

            for ($y=0; $y<20; $y++)
            {
                for ($x=0; $x<20; $x++)
                {
                    // $value = $x/20 * pi(); // for testing
                    $value = floatval($grid[$y][$x]) / $gridmax * pi();
                    $red   = intval(128 - cos($value + 0.5) * 127) * pow($value/pi(), 0.2);
                    $green = intval($value >= 1.5 ? (128 - cos(($value - 1.5)*2.3) * 127) : 0);
                    $blue  = intval(128 - cos(($value + 1.5)*1.8) * 127) * pow($value/pi(), 0.333);

                    echo "<span style=\"background-color: rgb($red, $green, $blue);\">&nbsp;&nbsp;&nbsp;&nbsp;</span>";
                }
                echo "<br>\n";
            }
            echo "</td>\n";
        }
        ?>
    </tr>
</table>

<?php

output_dlmenu_div();
