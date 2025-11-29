<?php
chdir(__DIR__);
require_once("../data/protutils.php");
require_once("../data/odorutils.php");
require_once("receptor_notes.php");
require_once("dlmenu.php");

$lrefs = [];

/*
foreach ($odors as $o)
{
    if (isset($o['activity'])) foreach ($o['activity'] as $a)
    {
        foreach ($a as $rcpid => $data)
        {
            if (@$data['adjusted_curve_top'] > 0) $prots[$rcpid]['has_agonists'] = true;
            else if (isset($data['type'])) switch ($data['type'])
            {
                case 'a': case 'vsa': case 'sa': case 'ma': case 'wa': case 'vwa':
                    $prots[$rcpid]['has_agonists'] = true;
                    break;
                case 'ia':
                    $prots[$rcpid]['has_inverse_agonists'] = true;
                    break;
                default:
                    ;
            }
        }
    }
} */

function get_refno($refurl)
{
    global $lrefs;
    foreach ($lrefs as $rno => $url) if ($url == $refurl) return $rno+1;
    $lrefs[] = $refurl;
    return count($lrefs);
}

$odor = find_odorant(@$_REQUEST['o']);
if (!$odor)
{
    header("Location: odorants.php");
    exit;
}

if (!file_exists($pqcorr = "../data/receptor_pq.json"))
{
    $correlations = correlate_receptors_aromanotes();

    $f = fopen($pqcorr, "wb");
    if (!$f) die("FAILED to open $pqcorr for writing.");
    fwrite($f, preg_replace("/([ \t]*)([^\\s]*) ([\\[{])\n/", "\$1\$2\n\$1\$3\n", json_encode($correlations, JSON_PRETTY_PRINT)));
    fclose($f);
}
else
{
    chdir(__DIR__);
    $correlations = json_decode(file_get_contents($pqcorr), true);
}

$predictions = [];
$predvals = [];
chdir(__DIR__);
$dock_results = json_decode(file_get_contents("../predict/dock_results.json"), true);
$odorname_under = str_replace(' ', '_', $odor["full_name"]);
foreach ($dock_results as $protid => $dr)
{
    if (isset($dr[$odorname_under]))
    {
        $abn = $ibn = 1;
        $abe = $ibe = 9.16e+58;
        if (isset($dr[$odorname_under]["ascores"]))
        {
            foreach ($dr[$odorname_under]["ascores"] as $aln => $ale) if ($ale < $abe)
            {
                $abn = $aln;
                $abe = $ale;
            }
        }
        if (isset($dr[$odorname_under]["iscores"]))
        {
            foreach ($dr[$odorname_under]["iscores"] as $iln => $ile) if ($ile < $ibe)
            {
                $ibn = $iln;
                $ibe = $ile;
            }
        }

        $predictions[$protid] = $dr[$odorname_under]["DockScore"];
        $predvals[$protid] =
        [
            "active_n" => $abn,
            "inactive_n" => $ibn,
            "Affinity" => $dr[$odorname_under]["Affinity"],
            "A100" => $dr[$odorname_under]["A100"],
            "Pose1" => $dr[$odorname_under]["a{$abn}_Pose1"],
            "LSP" => $dr[$odorname_under]["a{$abn}_LSP"],
            "LSRB" => $dr[$odorname_under]["a{$abn}_LSRB"],
            "Docker" => $dr[$odorname_under]["docker"],
            "Predate" => $dr[$odorname_under]["CalculateDate"]
        ];
    }
}

$page_title = $odor['full_name'];
$extra_js = ['js/tabs.js'];
$extra_css = ['assets/tabs.css'];

include("header.php");

?>
<style>
<?php output_dlmenu_style(); ?>
</style>
<script>
var viewer_loaded = false;
function load_viewer(obj)
{
    openTab(obj, 'Structure');
    if (!viewer_loaded)
    {
        window.setTimeout( function()
        {
            $('#viewer').on('load', function()
            {
                var embdd = $('#viewer')[0];
                $("[type=file]", embdd.contentDocument).hide();
                var filediv = $("#filediv", embdd.contentDocument)[0];
                $("#posey", embdd.contentDocument).css("position", "absolute").css("left", "-262144px");

                <?php
                $forms = [];
                if (@$odor['isomers']) $forms = array_merge($forms, $odor['isomers']);
                if (@$odor['forms']) $forms = array_merge($forms, $odor['forms']);

                if (count($forms))
                {
                    echo "filediv.innerHTML = \"".@$odor['full_name'];
                    if (!@$odor['isomers'])
                    {
                        echo " &#xb7; <small><a href=\\\"#\\\" onclick=\\\"load_remote_sdf('sdf.php?mol={$odor['oid']}');\\\">(normal)</a></small>";
                    }
                    foreach (array_keys($forms) as $form)
                    {
                        $formu = urlencode($form);
                        echo " &#xb7; <small><a href=\\\"#\\\" onclick=\\\"load_remote_sdf('sdf.php?mol={$odor['oid']}&form=$formu');\\\">$form</a></small>";
                    }
                    echo "\";";
                }
                else
                {
                    ?>filediv.innerText = "<?php echo @$odor['full_name']; ?>";
                <?php } ?>
            });
            $('#viewer')[0].src = '<?php echo "viewer.php?url=sdf.php&mol={$odor['oid']}"; ?>'; 
        }, 259); 
        viewer_loaded = true;
    }
}

window.setTimeout( function()
{
    $('#skeletal')[0].innerHTML = svg_from_smiles('<?php echo str_replace("\\", "\\\\", $odor['smiles']); ?>', 300, 300);
    $('#skeletal1')[0].innerHTML = svg_from_smiles('<?php echo str_replace("\\", "\\\\", $odor['smiles']); ?>', 300, 300);
    var boundary = parseInt($(".tab")[0].getClientRects()[0].bottom);
    $("#tabAroma").click();
}, 123);

<?php output_dlmenu_script(); ?>
</script>
<div class="tab" style="display: inline-block; margin-top: 30px;">
    <button class="tabstatic" id="tabFullName"><?php echo $odor['full_name']; ?></button>
	<button class="tablinks" id="tabAroma" onclick="openTab(this, 'Aroma');">Notes & Receptors</button>
    <?php if (count($predictions)) { ?>
	<button class="tablinks" id="tabPred" onclick="openTab(this, 'Predict');">Predictions</button>
    <?php } ?>
    <button class="tablinks" id="tabRefs" onclick="openTab(this, 'Refs');">References</button>
	<button	class="tablinks"
			id="tabStructure"
			onclick="load_viewer(this);"
			>3D Structure</button>
</div>

<div id="Aroma" class="tabcontent">

<div class="scrollw">
    <div>
        <div id="skeletal">&nbsp;</div>
        <img class="barchart" src="barchart.php?o=<?php echo urlencode($odor['smiles']); ?>">
    </div>
</div>

<p class="aromainfo">
    <?php
    for ($i=1; isset($odor["name$i"]); $i++)
    {
        if ($i > 1) echo ", ";
        echo $odor["name$i"];
    }
    if (isset($odor["iupac"]))
    {
        if ($i > 1) echo ", ";
        echo $odor["iupac"];
        if ($i == 1) echo "<br><br>";
    }
    if ($i > 1) echo "<br><br>";
    ?>
	<strong>SMILES:</strong><br>
	<?php echo $odor['smiles']; ?><br>
	<br>
    <strong>Aroma Description:</strong>
    <br>
    <?php 
    foreach ($odor['aroma'] as $refurl => $notes)
    {
        $comma = false;
        foreach ($notes as $note)
        {
            if ($comma) echo ", ";
            echo "$note";
            $comma = true;
        }
        $refno = get_refno($refurl);
        echo "<sup><a href=\"#\" onclick=\"openTab($('#tabRefs')[0], 'Refs');\">$refno</a></sup><br>";
    }

    if (@$odor['notes'])
    {
        echo "<br><strong>Notes:</strong><br>";
        echo "{$odor['notes']['text']}";

        $comma = false;
        echo "<sup>";
        foreach ($odor['notes']['refs'] as $refurl)
        {
            if ($comma) echo ", ";
            $refno = get_refno($refurl);
            echo "<a href=\"#\" onclick=\"openTab($('#tabRefs')[0], 'Refs');\">$refno</a>";
            $comma = true;
        }
        echo "</sup>";
        echo "<br>";
    }
    ?>
</p>

<div class="box" style="height: auto!important;">
<div class="row content scrollh">
<table class="rcplist">

<tr>
    <th>Receptor</th>
    <th>Expression</th>
    <th>log10 EC<sub>50</sub></th>
    <th>Adj. Top</th>
    <th>Antagonist?</th>
    <th>Correlated Perceptual Qualities</th>
</tr>

<?php

$sorted = [];
$tbltops = [];
$tblec50 = [];
$agonist = [];
$antagonist = [];
if (@$odor['activity']) foreach ($odor['activity'] as $refurl => $acv)
{
    $maxcurvtop = [];
    $minec50 = [];
    $refno = get_refno($refurl);
    foreach ($acv as $rcpid => $a)
    {
        $maxcurvtop[$rcpid] = false;
        $minec50[$rcpid] = false;

        if (@$a['antagonist']) $antagonist[$rcpid] = "Y";

        if (!isset($sorted[$rcpid])) $sorted[$rcpid] = 0.0;
        $ssamples = 0;
        if (!isset($a['adjusted_curve_top']) && @$a['type'] == 'na') $a['adjusted_curve_top'] = 0;
        else if (!isset($a['adjusted_curve_top']) && @$a['type'] == 'a') $a['adjusted_curve_top'] = "(agonist)";
        if (isset($a['adjusted_curve_top']))
        {
            if (!isset($tbltops[$rcpid])) $tbltops[$rcpid] = "";
            else $tbltops[$rcpid] .= ", ";

            $tbltops[$rcpid] .=
                is_numeric($a['adjusted_curve_top'])
                ? (round($a['adjusted_curve_top'], 4) . " <sup><a href=\"#\" onclick=\"openTab($('#tabRefs')[0], 'Refs');\">$refno</a></sup>")
                : ($a['adjusted_curve_top'] . " <sup><a href=\"#\" onclick=\"openTab($('#tabRefs')[0], 'Refs');\">$refno</a></sup><!-- \x50\x61\x79\x77\x61\x6c\x6c\x73\x20\x73\x75\x63\x6b\x2e -->");
            $sorted[$rcpid] += is_numeric($a['adjusted_curve_top']) ? $a['adjusted_curve_top'] : 10;
            $ssamples++;

            if (false===$maxcurvtop[$rcpid] || $a['adjusted_curve_top'] > $maxcurvtop[$rcpid]) $maxcurvtop[$rcpid] = $a['adjusted_curve_top'];
        }
        else if (isset($a['type']))
        {
            $top = false; $tosort = 0;
            if ($a['type'] == "na") $top = "0";
            else if ($a['type'] == "vsa" || $a['type'] == "sa" || $a['type'] == "ma" || $a['type'] == "wa" || $a['type'] == "vwa")
            {
                $top = "(agonist)";
                if ($a['type'] == "vsa") $tosort = 10;
                if ($a['type'] == "sa") $tosort = 8;
                if ($a['type'] == "ma") $tosort = 5;
                if ($a['type'] == "wa") $tosort = 1;
                if ($a['type'] == "vwa") $tosort = 0.1;
            }
            else if ($a['type'] == "pa")
            {
                $top = "(probable&nbsp;agonist)";
                $tosort = 0.001;
            }

            if ($top)
            {
                if (!isset($tbltops[$rcpid])) $tbltops[$rcpid] = "";
                else $tbltops[$rcpid] .= ", ";
                $tbltops[$rcpid] .= $top . " <sup><a href=\"#\" onclick=\"openTab($('#tabRefs')[0], 'Refs');\">$refno</a></sup>";
                $sorted[$rcpid] += $tosort;
                $ssamples++;
            }
        }
        if (isset($a['ec50']))
        {
            if (!isset($tblec50[$rcpid])) $tblec50[$rcpid] = "";
            else $tblec50[$rcpid] .= ", ";

            $tblec50[$rcpid] .= round($a['ec50'], 4) . " <sup><a href=\"#\" onclick=\"openTab($('#tabRefs')[0], 'Refs');\">$refno</a></sup>";
            if (!isset($a['adjusted_curve_top']) || floatval($a['adjusted_curve_top']) > 0)
            {
              $sorted[$rcpid] -= $a['ec50']*2;
              $ssamples++;
            }

            if (false===$minec50[$rcpid] || $a['ec50'] < $minec50[$rcpid]) $minec50[$rcpid] = $a['ec50'];
        }
        if ($ssamples) $sorted[$rcpid] /= $ssamples;
    }

    foreach (array_keys($maxcurvtop) as $rcpid)
    {
        if ($maxcurvtop[$rcpid] > 0) $agonist[$rcpid] = true;
        if ($maxcurvtop[$rcpid] >= 0 && $minec50[$rcpid] < 0) $agonist[$rcpid] = true;
    }
}

arsort($sorted);

$isagonist = [];
foreach (array_keys($sorted) as $rcpid)
{
    echo "<tr>\n";
    echo "<td><a href=\"receptor.php?r=$rcpid\">$rcpid</a></td>\n";
    echo "<td>".(@$prots[$rcpid]['expression'] ?: '?')."</td>\n";

    echo "<td style=\"white-space: nowrap;\">" . $dispec50 = (@$tblec50[$rcpid] ?: "-") . "</td>\n";
    echo "<td style=\"white-space: nowrap;\">" . $disptop =
    (
        (floatval(@$tbltops[$rcpid]) || !floatval(@$tblec50[$rcpid]))
        ? @$tbltops[$rcpid]
        : "-"
    ) . "</td>\n";

    if (@$antagonist[$rcpid]) echo "<td>Y</td>";
    else
    {
        echo "<td>&nbsp;</td>";
        $isagonist[$rcpid] = (floatval(@$tbltops[$rcpid]) > 0 || floatval(@$tblec50[$rcpid]) < 0);
    }

    if (@$agonist[$rcpid])
    {
        $notes = substr(get_notes_for_receptor($rcpid, $correlations), 0, 123);
        $notes = make_clickable_notes(array_slice(explode(", ", $notes), 0, 10));
        $notes = implode(", ", $notes);
        if (substr($notes, 0, 1) == '(')
        {
            $notes = "<i class=\"dim\">$notes</i>";
        }
        echo "<td style=\"white-space: nowrap;\">$notes</td>\n";
    }
    else
    {
        echo "<td>&nbsp;</td>\n";
    }    

    echo "</tr>\n";
}

?>
</table>
</div>
</div>

</div>
<div id="Predict" class="tabcontent">

<div class="scrollw">
    <div>
        <div id="skeletal1">&nbsp;</div>
        <img class="barchart" src="barchart.php?m=p&o=<?php echo urlencode($odor['smiles']); ?>&t=<?php echo time(); ?>">
    </div>
</div>

<p class="aromainfo">
    <?php
    for ($i=1; isset($odor["name$i"]); $i++)
    {
        if ($i > 1) echo ", ";
        echo $odor["name$i"];
    }
    if (isset($odor["iupac"]))
    {
        if ($i > 1) echo ", ";
        echo $odor["iupac"];
        if ($i == 1) echo "<br><br>";
    }
    if ($i > 1) echo "<br><br>";
    ?>
	<strong>SMILES:</strong><br>
	<?php echo $odor['smiles']; ?><br>
	<br>
    <strong>Aroma Description:</strong>
    <br>
    <?php 
    foreach ($odor['aroma'] as $refurl => $notes)
    {
        $comma = false;
        foreach ($notes as $note)
        {
            if ($comma) echo ", ";
            echo "$note";
            $comma = true;
        }
    }

    if (@$odor['comment']) echo "<br>{$odor['comment']}<br>";
    ?>
</p>

<div class="box" style="height: auto!important;">
<div class="row content scrollh50">
<table class="rcplist">

<tr>
    <th>Receptor</th>
    <th width="5%">Expr.%</th>
    <th width="5%">Agonist?</th>
    <th>Dock Score</th>
    <th>Known agonist</th>
    <th>Correlated Perceptual Qualities</th>
</tr>

<?php

arsort($predictions);

foreach ($predictions as $rcpid => $score)
{
    $aep = all_empirical_pairs_for_receptor($rcpid, true, true);
    $orphan = count($aep) ? false : true;
    echo "<tr>\n";
    $astyle = $orphan ? "font-style: italic;" : "";
    echo "<td><a href=\"receptor.php?r=$rcpid\" style=\"$astyle\">$rcpid</a></td>\n";
    echo "<td>".(@$prots[$rcpid]['expression'] ?: '?')."</td>\n";

    $aff = $predvals[$rcpid]["Affinity"];
    $Pose1 = $predvals[$rcpid]["Pose1"];
    $LSP = round($predvals[$rcpid]["LSP"], 4);
    $LSRB = round($predvals[$rcpid]["LSRB"], 4);
    $LSdelta = $LSRB - $LSP;
    $abn = $predvals[$rcpid]["active_n"];
    $ibn = $predvals[$rcpid]["inactive_n"];

    echo "<td style=\"white-space: nowrap;\">".(@$isagonist[$rcpid] ? "&#x2714;" : "&nbsp;")."</td>\n";
    $lig = urlencode(str_replace(' ', '_', $odor["full_name"]));
    echo "<td><a href=\"#\" onclick=\"show_dlmenu(event, '$rcpid', '$lig', '{$predvals[$rcpid]['Predate']}', '{$predvals[$rcpid]['Docker']}', '$abn', '$ibn');\">"
        .round($score, 4)."</a></td>";
    // echo "<td style=\"white-space: nowrap;\">$score</td>\n";

    if ($score > 0 && count($aep))
    {
        $agid = array_keys($aep)[0];
        $ag = $odors[$agid];
        $agname = str_replace(' ', '&nbsp;', $ag["full_name"]);
        echo "<td style=\"white-space: nowrap;\">$agname&nbsp;&nbsp;&nbsp;</td>\n";
    }
    else
    {
        echo "<td>&nbsp;</td>\n";
    }

    if ($score > 0)
    {
        $notes = substr(get_notes_for_receptor($rcpid, $correlations), 0, 203);
        $notes = make_clickable_notes(explode(", ", $notes));
        $notes = implode(", ", $notes);
        if (substr($notes, 0, 1) == '(')
        {
            $notes = "<i class=\"dim\">$notes</i>";
        }
        echo "<td style=\"white-space: nowrap;\">$notes</td>\n";
    }
    else
    {
        echo "<td>&nbsp;</td>\n";
    }

    echo "</tr>\n";
}

?>
</table>
</div>
</div>

<p class="small">
Dock Score is a measure of how strongly the algorithm thinks the odorant is likely to be an agonist of the receptor.<br>
Receptors in italics are "orphans", i.e. receptors whose agonists have not been identified experimentally.
</div>
</p>

<div id="Refs" class="tabcontent">
<div class="box">
<div class="row content scrollh">
<?php
foreach ($lrefs as $idx => $refurl)
{
    echo "<a href=\"$refurl\"><p>\n";
    $idx1 = $idx + 1;
    echo "$idx1.) ";
    if (!@$refs[$refurl]['citation'] && false!==strpos($refurl, "leffingwell.com")) $refs[$refurl]['citation'] = "Leffingwell";
    echo @$refs[$refurl]['citation'] ?: $refurl;
    echo "</p></a>\n";
}
?>
<!-- p>References for aroma perceptual qualities should not be taken to indicate that the authors of outside studies necessarily
assigned aroma notes to the neurons that receive input from any given receptor.
Rather, the findings of outside studies often constitute the information on which we base our own perceptual quality assignments.
</p -->
</div>
</div>
</div>

<div id="Structure" class="tabcontent">
    <iframe id="viewer"></iframe>
</div>

<?php output_dlmenu_div();