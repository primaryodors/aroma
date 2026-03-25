<?php
chdir(__DIR__);
require_once("../data/protutils.php");
require_once("../data/odorutils.php");
require_once("../data/statistics.php");
require_once("dlmenu.php");

$extra_js = ['js/tabs.js'];
$extra_css = ['assets/tabs.css'];
$page_title = "AROMA Docks";
include("header.php");

chdir(__DIR__);

$celph = $_SERVER['SERVER_NAME'];

$fh = 44;
$fw = 30;
$nih = 25;
$nt = 20;
$noh = $nih + $nt;
$nw = 18;
$filter_svgdat = "m 0,0 $fw,$fh 0,$nih $nw,$nt 0,-$noh $fw,-$fh Z";


$dockbsr = [];
$dbsrfn = "../out/dockbsr.json";
if (file_exists($dbsrfn)) $dockbsr = json_decode(file_get_contents($dbsrfn), true);

?>
<style>
<?php output_dlmenu_style(); ?>
.hop
{
    text-align: center!important;
}
.modeL
{
    background-color: rgba(0, 128, 64, 0.08);
    text-align: center!important;
}
.modeS2
{
    background-color: rgba(128, 128, 0, 0.08);
    text-align: center!important;
}
.modeI2
{
    background-color: rgba(128, 40, 0, 0.08);
    text-align: center!important;
}
.modeO1
{
    background-color: rgba(128, 0, 96, 0.08);
    text-align: center!important;
}
.modeactive
{
    background-color: rgba(64, 128, 64, 0.08);
    text-align: center!important;
}
.modeinactive
{
    background-color: rgba(0, 64, 128, 0.08);
    text-align: center!important;
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

<?php
$allmodes = [];
$protfs = [];
foreach ($prots as $protid => $p)
{
    // if (isset($_REQUEST['r']) && $protid != $_REQUEST['r']) continue;
    $fam = family_from_protid($protid);
    $dockpath = "../out/$fam/$protid";
    if (!file_exists($dockpath)) continue;
    $dir = dir($dockpath);
    $files = [];
    while (false!==($fname=$dir->read())) $files[] = $fname;
    natsort($files);
    $protfs[$protid] = $files;

    foreach ($protfs[$protid] as $fname)
    {
        if (substr($fname, -5) != ".dock") continue;
        if (false===strpos($fname, "~")) continue;
        list($odor, $mode, $opfisehciet) = explode('.', explode('~', $fname)[1], 3);
        if ($mode == "L") $mode = "000000000";
        else if ($mode == "S2") $mode = "000000001";
        $allmodes[$mode] = $mode;
    }
}
$allmodes = array_values($allmodes);
natsort($allmodes);
foreach ($allmodes as $k => $v)
{
    if ($v == "000000000") $allmodes[$k] = "L";
    if ($v == "000000001") $allmodes[$k] = "S2";
}
?>

<table class="liglist">
    <tr><th>Receptor</th>
        <th>Odorant</th>
        <th>Dock Date</th>
        <th colspan="<?php echo count($allmodes); ?>" class="hop">Enthalpy, Occlusion, # Poses</th>
        <th>Agonist Y/N,</td>
        <th>Prediction</td>
    </tr>
    <tr><th>ID</th>
        <th>Name</th>
        <th>(UTC)</th>
        <!-- th>Dock Energies</th>
        <th>Occlusion</th>
        <th>Poses</th -->
        <?php
        foreach ($allmodes as $m)
        {
            if (strlen($m) > 3)
                echo "<th class=\"mode$m\">$m</th>\n";
            else
                echo "<th class=\"mode$m\">+GNA$m</th>\n";
        }
        ?>
        <th>EC<sub>50</sub>, crvtop</td>
        <th>Score</td>
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

$poke = (false!==strpos($prod, "\x70\x72\x69\x6d\x61\x72\x79"))
    || (false!==strpos($prod, "\x75\x6d\x6f\x70"))
    || (false!==strpos($prod, "\x2e\x6f\x72\x67"))
    || (false!==strpos($prod, "\x2e\x6e\x65\x74"));

$graphdat =
[
    0 => [],
    1 => [],
    2 => [],
    3 => []
];
$right = $wrong = 0;
$topxy = $ec50xy = ["x" => [], "y" => []];
foreach ($prots as $protid => $p)
{
    if (isset($_REQUEST['r']) && $protid != $_REQUEST['r']) continue;
    $fam = family_from_protid($protid);

    $rows = [];
    $modes = [];
    $dates = [];
    $incplt = [];
    foreach ($protfs[$protid] as $fname)
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
        if (!isset($modes[$rowid])) $modes[$rowid] = [];
        $modes[$rowid][] = $mode;

        $fpn = "../out/$fam/$protid/$fname";
        $fmt = filemtime($fpn);

        if (!isset($dates[$rowid])) $dates[$rowid] = $fmt;
        else $dates[$rowid] = max($dates[$rowid], $fmt);

        set_time_limit(600);
        if (isset($cached[$fname]) /*&& intval($cached[$fname]['nump'])*/ && $fmt < $cachemt)
        {
            extract($cached[$fname]);
            // echo "<p>$protid $odor $mode<br><pre>".print_r($cached[$fname], true)."</pre></p>";
            ob_flush(); flush();
        }
        else
        {
            ob_flush(); flush();
            echo "<!-- Reading $fpn ... -->\n";
            chdir(__DIR__);
            $c = file_get_contents($fpn);
            if (strlen($c) < 500)
            {
                $incplt[$rowid] = true;
                continue;
            }
            $lines = explode("\n", $c);
            $benerg = 0;
            $lsfe = $lsbe = $phfe = $phbe = $ewde = 0;
            $nump = 0;
            $occl = 0;
            $tds = 0;
            foreach ($lines as $ln) 
            {
                if (substr($ln, 0, 21) == "Disqualified because:" || substr($ln, 0, 25) == "Disqualified for reasons.")
                {
                    $benerg = 1000;
                    $lsfe = $lsbe = $phfe = $phbe = $ewde = $occl = $tds = $nump = 0;
                    break;
                }
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
                else if (!$ewde && substr($ln, 0, 37) == "Estimated water displacement energy: ")
                {
                    $ewde = floatval(substr($ln, 37));
                }
                else if (!$occl && substr($ln, 0, 25) == "Ligand pocket occlusion: ")
                {
                    $occl = floatval(substr($ln, 25));
                }
                else if (!$tds && substr($ln, 0, 19) == "Estimated T_Delta_S: ")
                {
                    $tds = floatval(substr($ln, 19));
                }
                else if (substr($ln, 0, 6) == "Pose: ") $nump++;

                if ($mode == "active" && $nump < 2 && false!==strpos($ln, "~(ligand)"))
                {
                    list($lcntct, $strength) = explode(": ", $ln, 2);
                    $strength = floatval($strength);
                    if ($strength <= -0.5)
                    {
                        $liga = false;
                        $resno = false;
                        list($I, $II) = explode('~', $lcntct);
                        list($res, $resa) = explode(':', $I);
                        list($lig, $liga) = explode(':', $II);
                        if (!$liga) continue;
                        $resno = intval(preg_replace("/[^0-9]/", "", $res));
                        if (!$resno) continue;
                        $bw = bw_from_resno($protid, $resno);
                        if (!isset($dockbsr[$protid])) $dockbsr[$protid] = [];
                        if (!isset($dockbsr[$protid][$odor])) $dockbsr[$protid][$odor] = [];
                        $dockbsr[$protid][$odor][$bw] = $liga;
                    }
                }
            }
            if (!$nump) $nump = "-";

            if (0)              // enable this and delete docklist.cache.json to update occlusions whenever changed the occlusion calculation.
            {
                chdir(__DIR__);
                $fpnrel = substr($fpn, 1);
                $cmd = "../test/occlusion_test \"$fpnrel\"";
                $result = [];
                exec($cmd, $result);
                foreach ($result as $ln)
                {
                    $pieces = explode(": ", $ln);
                    if ($pieces[0] == "Occlusion") $occl = floatval(@$pieces[1]);
                }
            }

            $cached[$fname] =
            [
                "benerg" => $benerg,
                "lsfe" => $lsfe,
                "lsbe" => $lsbe,
                "phfe" => $phfe,
                "phbe" => $phbe,
                "ewde" => $ewde,
                "occl" => $occl,
                "nump" => $nump,
                "tds"  => $tds,
            ];
        }

        $rows[$rowid]["benerg_raw_$mode"] = $benerg;
        $rows[$rowid]["benerg_$mode"] = $benerg + ($lsbe - $lsfe) + ($phbe - $phfe) + $ewde - $tds;
        $rows[$rowid]["nump_$mode"] = $nump;
        $rows[$rowid]["occl_$mode"] = $occl;

        $agonist = "?";
        $pair = best_empirical_pair($protid, $odor, true);
        if ($pair)
        {
            if (isset($pair["adjusted_curve_top"])) $rows[$rowid]["top"] = round($pair["adjusted_curve_top"], 3);
            if (isset($pair["ec50"])) $rows[$rowid]["ec50"] = $pair["ec50"];
            if (isset($pair["adjusted_curve_top"]))
            {
                $act = floatval($pair["adjusted_curve_top"]);
                if ($act > 0) $agonist = "Y";
                else if (!$act) $agonist = "N";
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

                    case "pa":
                        $agonist = "Y?";
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

    foreach ($rows as $k => $r)
    {
        list($protid, $odor) = explode("~", $k);
        $benerg_active = 0;
        // print_r($r);
        $benerg_ = []; $nump_ = []; $occl_ = []; $benerg_raw_ = [];
        foreach ($modes[$k] as $mode)
        {
            $benerg_[$mode] = @$r["benerg_$mode"] ?: 0;
            $nump_[$mode] = @$r["nump_$mode"] ?: 0;
            $occl_[$mode] = @$r["occl_$mode"] ?: 0;
            $benerg_raw_[$mode] = @$r["benerg_raw_$mode"] ?: 0;
            if ($mode != "inactive" && (!$benerg_active || $benerg_active > $benerg_[$mode]))
            {
                $benerg_active = $benerg_[$mode];
                $nump_active = $nump_[$mode];
                $occl_active = $occl_[$mode];
                $benerg_raw_active = $benerg_raw_[$mode];
            }
        }
        $top = $ec50 = false;
        extract($r);

        echo "<tr>\n";

        echo "<td><a href=\"receptor.php?r=$protid\">$protid</a>";
        if ($frcp != $protid)
            echo " <a href=\"docklist.php?r=$protid\"><svg height=\"13px\" viewBox=\"0 0 80 90\" xmlns=\"http://www.w3.org/2000/svg\"><path fill=\"#50cea8\" d=\"$filter_svgdat\"></path></svg></a>";
        echo "</td>\n";
        $frcp = $protid;

        $o = find_odorant($odor);
        $fn = $o['full_name'];
        $fnu = urlencode(str_replace(' ', '_', $fn));
        echo "<td><a href=\"odorant.php?o={$o['oid']}\">$fn</a>";
        if ($flig != $o['oid'])
            echo " <a href=\"docklist.php?o={$o['oid']}\"><svg height=\"13px\" viewBox=\"0 0 80 90\" xmlns=\"http://www.w3.org/2000/svg\"><path fill=\"#50cea8\" d=\"$filter_svgdat\"></path></svg></a>";
        echo "</td>\n";
        $flig = $o['oid'];

        echo "<td>".date("Y-m-d H:i:s", $dates[$k])."</td>";

        foreach ($allmodes as $m)
        {
            echo "<td class=\"mode$m\">";
            if (@$benerg_[$m]) echo "<a href=\"viewer.php?view=dock&prot=$protid&odor=$fnu&mode=$m\" target=\"_dock\">";
            if (@$benerg_[$m] >= 200) echo "(fail)";
            else
            {
                echo @$benerg_[$m] ? round(@$benerg_[$m], 2) : "-";
                if (@$benerg_[$m])
                {
                    echo " / ";
                    echo @$occl_[$m] ? round($occl_[$m], 2) : "-";
                    echo " / ";
                    echo @$nump_[$m] ?: "-";
                }
            }
            if (@$benerg_[$m]) echo "</a>";
            echo "</td>";
        }

        /*
        echo "<td>";

        if ($benerg_active)
        {
            if (@$benerg_["active"])
            {
                if ($benerg_active >= 200) $dispe = "(fail)";
                else $dispe = round($benerg_active, 4);
                echo "<a href=\"viewer.php?view=dock&prot=$protid&odor=$fnu&mode=active\" target=\"_dock\">$dispe</a>";
            }
            else if (count($modes[$k]))
            {
                $frist = true;
                foreach ($modes[$k] as $mode)
                {
                    if ($mode == "inactive") continue;
                    if (!$frist) echo ", ";
                    if (strlen($mode) < 6) echo "$mode:";
                    if ($benerg_[$mode] >= 200) $dispe = "(fail)";
                    else $dispe = round($benerg_[$mode], 4);
                    if ($poke) echo "<span class=\"color: #f00;\">$dispe</span>";
                    else echo "<a href=\"viewer.php?view=dock&prot=$protid&odor=$fnu&mode=$mode\" target=\"_dock\">$dispe</a>";
                    $frist = false;
                }
            }
            else echo "-";
        }
        else echo "-";

        echo " / ";
        if ($benerg_inactive)
        {
            if ($benerg_inactive >= 200) $dispe = "(fail)";
            else $dispe = round($benerg_inactive, 4);
            echo "<a href=\"viewer.php?view=dock&prot=$protid&odor=$fnu&mode=inactive\" target=\"_dock\">$dispe</a>";
        }
        else echo "-";
        echo "</td>";

        echo @"<td>";

        if ($occl_["active"])
            echo round($occl_active, 3) ?: "-";
        else
        {
            $frist = true;
            foreach ($modes[$k] as $mode)
            {
                if ($mode == "inactive") continue;
                if (!$frist) echo ", ";
                if (strlen($mode) < 6) echo "$mode:";
                echo $occl_[$mode] ? round($occl_[$mode], 3) : "-";
                $frist = false;
            }
        }

        echo " / " . (round($occl_inactive, 3) ?: "-") . "</td>\n";


        echo @"<td>" . ($nump_active ?: "-") . " / " . ($nump_inactive ?: "-") . "</td>\n";
        */

        $dtop = $top ?: "-";
        $dec50 = $ec50 ?: "-";
        echo @"<td>$agonist $dec50/$dtop</td>\n";

        if ($benerg_inactive > 0) $benerg_inactive = 0;
        if ($occl_active >= 0.6)
            $prediction = max(0, -$benerg_raw_active) * 3.0
                * equilibrium(-$benerg_active * $occl_active, $nump_inactive ? (-$benerg_inactive * $occl_inactive) : 0);
        else $prediction = 0;

        if ($top)
        {
            $topxy['x'][] = floatval($top);
            $topxy['y'][] = $prediction;
        }
        if ($ec50)
        {
            $ec50xy['x'][] = floatval($ec50);
            $ec50xy['y'][] = $prediction;
        }

        $prediction = round($prediction, 2);

        $color = "";
        if (@$incplt[$rowid])
        {
            $prediction = "...";
        }
        else
        {
            if ($agonist == 'Y')
            {
                $graphdat[0][] = [$benerg_active, $occl_active, $prediction];
                $graphdat[2][] = [$benerg_active - $benerg_inactive, $occl_active, $prediction];
                if ($prediction > 0)
                {
                    $color = "color: #0c0;";
                    $right++;
                }
                else if ($top && $top < 2.5)
                {
                    $color = "color: #f90;";
                }
                else if (!$top && $ec50 > -3)
                {
                    $color = "color: #f90;";
                }
                else
                {
                    $color = "color: #f00;";
                    $wrong++;
                }
            }
            if ($agonist == 'Y?')
            {
                if ($prediction > 0)
                {
                    $color = "color: #0cc;";
                    $right++;
                }
                else
                {
                    $color = "color: #960;";
                    $wrong++;
                }
            }
            else if ($agonist == 'inv' || $agonist == 'ant')
            {
                $graphdat[1][] = [$benerg_active, $occl_active, $prediction];
                $graphdat[3][] = [$benerg_active - $benerg_inactive, $occl_active, $prediction];
                if ($prediction <= 0)
                {
                    $color = "color: #0c0;";
                    $right++;
                }
                else if ($prediction < 4)
                {
                    $color = "color: #f90;";
                }
                else
                {
                    $color = "color: #f00;";
                    $wrong++;
                }
            }
            else if ($agonist == 'N')
            {
                $graphdat[1][] = [$benerg_active, $occl_active, $prediction];
                $graphdat[3][] = [$benerg_active - $benerg_inactive, $occl_active, $prediction];
                if (!$prediction)
                {
                    $color = "color: #0c0;";
                    $right++;
                }
                else if ($prediction < 8)
                {
                    $color = "color: #f90;";
                }
                else
                {
                    $color = "color: #f00;";
                    $wrong++;
                }
            }
        }

        echo "<td style=\"$color\">$prediction</td>\n";

        $benerg_active = $occl_active = $nump_active = $benerg_inactive = $occl_inactive = $nump_inactive = $prediction = 0;
    }
}

chdir(__DIR__);
file_put_contents($cachefn, json_encode_pretty($cached));
file_put_contents($dbsrfn, json_encode_pretty($dockbsr));

?></table>
<?php 
if ($right + $wrong)
{
    echo "<p>Accuracy: " . round(100.0 * floatval($right) / ($right+$wrong), 2) . "%<br>";
    echo "Correlations:";
    if (count($topxy['x']) > 3) echo " Top = " . round(correlationCoefficient($topxy['x'], $topxy['y']), 4);
    if (count($ec50xy['x']) > 3) echo " EC<sub>50</sub> = " . round(correlationCoefficient($ec50xy['x'], $ec50xy['y']), 4);
    echo "</p>"; 
}
?>

<?php if (0) { ?>
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

<?php }

output_dlmenu_div();
