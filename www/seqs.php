<?php

chdir(__DIR__);
require_once("../data/protutils.php");
require_once("../data/odorutils.php");

if (file_exists("tags.php")) include("tags.php");

include("header.php");

chdir(__DIR__);
chdir("..");

$seqlns = explode("\n", file_get_contents("data/sequences_aligned.txt"));
$seqa = [];
$ells = [];
$emms = [];
$rgns = [];
$rgne = [];
$rg50 = [];
$lells = [];
$lemms = [];
$lrgns = [];
$lrgne = [];
$lrg50 = [];
$big5 = ["OR1A1", "OR1D2", "OR1G1", "OR2W1", "OR52D1"];
foreach ($seqlns as $ln)
{
    $c = substr($ln, 0, 1);
    if ($c >= 'A' && $c <= 'Z')
    {
        $protid = explode(' ',$ln)[0];
        $seqa[$protid] = $ln;
        $fam = family_from_protid($protid);
        if (count($lells) && !isset($ells[$fam])) $ells[$fam] = $lells;
        if (count($lemms) && !isset($emms[$fam])) $emms[$fam] = $lemms;
        if (count($lrg50) && !isset($rg50[$fam])) $rg50[$fam] = $lrg50;
        if (count($lrgne) && !isset($rgne[$fam])) $rgne[$fam] = $lrgne;
        if (count($lrgns) && !isset($rgns[$fam])) $rgns[$fam] = $lrgns;
    }
    else if (preg_match("/^[LM ]+$/", $ln))
    {
        $lells = [];
        $lemms = [];
        $i=0;
        while ($i = strpos($ln, 'L', $i+1)) $lells[$i] = true;
        $i=0;
        while ($i = strpos($ln, 'M', $i+1)) $lemms[$i] = true;
    }
    else if (preg_match("/^[A-Z0-9\\| -]+$/", $ln))
    {
        $ln .= " ";
        $lrgns = [];
        $lrgne = [];
        $lrg50 = [];
        $i=0;
        while ($i = strpos($ln, '|', $i+1)) $lrg50[$i] = true;
        $l = strlen($ln);
        $currrgn = "";
        for ($i=0; $i<$l; $i++)
        {
            $c = substr($ln, $i, 1);
            if ($c != ' ')
            {
                if (!$currrgn && $c != "|")
                {
                    $currrgn = $c;
                    $lrgns[503] = $i;
                }
                else
                {
                    if ($c >= 'A' && $c <= 'Z') $currrgn .= $c;
                    else if (is_numeric($c)) $currrgn .= $c;
                }
            }
            else if ($currrgn)
            {
                $lrgne[$currrgn] = $i-1;
                $lrgns[$currrgn] = $lrgns[503];
                unset($lrgns[503]);
                $currrgn = "";
            }
        }
    }
}

foreach ($odors as $o)
{
    if (isset($o['activity'])) foreach ($o['activity'] as $a)
    {
        foreach ($a as $rcpid => $data)
        {
            if (@$data['adjusted_curve_top'] > 0) $prots[$rcpid]['has_agonists'] = true;
            else if (isset($data['type'])) switch ($data['type'])
            {
                case 'vsa': case 'sa': case 'ma': case 'wa': case 'vwa':
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
}

$tree = [];
foreach ($prots as $protid => $p)
{
    if (!isset($p['btree'])) continue;
    $path = $p['btree'];
    if (!isset($tree[$path])) $tree[$path] = [];
    $tree[$path] = $protid;
}

ksort($tree, SORT_STRING);

?>
<style>
#seqdiv
{
    width: 100%;
    display: block;
    height: 100%; 
    overflow: scroll;
}

#seqdiv pre
{
    padding: 10px;
    background-color: #11161c;
    color: #89a;
    width: fit-content;
}

#seqdiv .big5
{
    color: #ccc;
}

#seqdiv .rg50
{
    color: #fff;
    background-color: #600;
}

#seqdiv .bsr
{
    font-weight: bold;
    color: #9df;
    background-color: #107;
}

#seqdiv .mcoord
{
    font-weight: bold;
    color: #600;
    background-color: #fd6;
}

#seqdiv .rgn1
{
    color: #bbb;
}

#seqdiv .rgn2
{
    color: #9bd;
}

#seqdiv .rgn3
{
    color: #7ba;
}

#seqdiv .rgn4
{
    color: #ac9;
}

#seqdiv .rgn5
{
    color: #bb8;
}

#seqdiv .rgn6
{
    color: #b97;
}

#seqdiv .rgn7
{
    color: #d89;
}

#seqdiv .rgn8
{
    color: #98c;
}

.OR1
{
    background-color: #1e2424;
}

.OR2
{
    background-color: #21221c;
}

.OR3
{
    background-color: #17281c;
}

.OR4
{
    background-color: #1a1633;
}

.OR5
{
    background-color: #261a23;
}

.OR6
{
    background-color: #11222a;
}

.OR7
{
    background-color: #26161c;
}

.OR8
{
    background-color: #281e1c;
}

.OR9
{
    background-color: #291a1c;
}

.OR10
{
    background-color: #262022;
}

.OR11
{
    background-color: #1b212c;
}

.OR12
{
    background-color: #1d281c;
}

.OR13
{
    background-color: #112425;
}

.OR14
{
    background-color: #112033;
}

.OR51
{
    background-color: #202321;
}

.OR52
{
    background-color: #21162c;
}

.OR56
{
    background-color: #11261c;
}

.TAAR
{
    background-color: #15163c;
}

.VN1R
{
    background-color: #21262c;
}

</style>
<h1>Sequence Alignments</h1>

<div id="seqdiv">
    <pre><?php // print_r($rgns); exit;
            $ffam = "";
            foreach ($tree as $protid)
            {
                $fam = family_from_protid($protid);
                if ($ffam != $fam)
                {
                    $lells = @$ells[$fam];
                    $lemms = @$emms[$fam];
                    $lrg50 = @$rg50[$fam];
                    $lrgns = array_flip(@$rgns[$fam] ?: []);
                    $lrgne = array_flip(@$rgne[$fam] ?: []);
                    if ($ffam) echo "\n";
                }

                $lseqln = @$seqa[$protid];

                $ldispln = "";
                $nemms = 0;
                $currrgn = "";
                for ($i=0; $i<strlen($lseqln); $i++)
                {
                    if (@$lrgns[$i])
                    {
                        $currrgn = $lrgns[$i];
                        $rgno = intval(substr($currrgn, 3));
                        $ldispln .= "<span class=\"rgn$rgno\">";
                    }
                    else if (@$lrgne[$i])
                    {
                        $currrgn = "";
                        $ldispln .= "</span>";
                    }

                    $c = substr($lseqln, $i, 1);

                    if (@$lrg50[$i]) $ldispln .= "<span class=\"rg50\">$c</span>";
                    else if (@$lemms[$i] && false!==strpos("CHM", $c))
                    {
                        $ldispln .= "<span class=\"mcoord\">$c</span>";
                        $nemms++;
                    }
                    else if (@$lells[$i] || @$lemms[$i]) $ldispln .= "<span class=\"bsr\">$c</span>";
                    else if ($c == ' ') $ldispln .= '.';
                    else $ldispln .= $c;
                }
                if ($currrgn) $ldispln .= "</span>";

                if ($nemms < 3) $ldispln = str_replace("class=\"mcoord\"", "class=\"bsr\"", $ldispln);

                $ex = floatval(@$prots[$protid]["expression"]);
                if ($ex >= 98) $ldispln = str_replace("$protid...", "$protid**.", $ldispln);
                else if ($ex >= 75) $ldispln = str_replace("$protid...", "$protid*..", $ldispln);

                if (!@$prots[$protid]['has_agonists']) $ldispln = str_replace($protid, "<i>$protid</i>", $ldispln);
                else if (in_array($protid, $big5)) $ldispln = str_replace($protid, "<span class=\"big5\">$protid</span>", $ldispln);

                echo "<span class=\"$fam\">$ldispln</span>\n";
                $ffam = $fam;
            }
        ?>

Legend:
<span class="OR1">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR2">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR3">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR4">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR5">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR6">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR7">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR8">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR9">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR10">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR11">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR12">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR13">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR14">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR51">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR52">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR56">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="TAAR">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="VN1R">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span>
<span class="OR1">&nbsp;&nbsp;&nbsp;OR1&nbsp;&nbsp;&nbsp;</span> <span class="OR2">&nbsp;&nbsp;&nbsp;OR2&nbsp;&nbsp;&nbsp;</span> <span class="OR3">&nbsp;&nbsp;&nbsp;OR3&nbsp;&nbsp;&nbsp;</span> <span class="OR4">&nbsp;&nbsp;&nbsp;OR4&nbsp;&nbsp;&nbsp;</span> <span class="OR5">&nbsp;&nbsp;&nbsp;OR5&nbsp;&nbsp;&nbsp;</span> <span class="OR6">&nbsp;&nbsp;&nbsp;OR6&nbsp;&nbsp;&nbsp;</span> <span class="OR7">&nbsp;&nbsp;&nbsp;OR7&nbsp;&nbsp;&nbsp;</span> <span class="OR8">&nbsp;&nbsp;&nbsp;OR8&nbsp;&nbsp;&nbsp;</span> <span class="OR9">&nbsp;&nbsp;&nbsp;OR9&nbsp;&nbsp;&nbsp;</span> <span class="OR10">&nbsp;&nbsp;&nbsp;OR10&nbsp;&nbsp;&nbsp;</span> <span class="OR11">&nbsp;&nbsp;&nbsp;OR11&nbsp;&nbsp;&nbsp;</span> <span class="OR12">&nbsp;&nbsp;&nbsp;OR12&nbsp;&nbsp;&nbsp;</span> <span class="OR13">&nbsp;&nbsp;&nbsp;OR13&nbsp;&nbsp;&nbsp;</span> <span class="OR14">&nbsp;&nbsp;&nbsp;OR14&nbsp;&nbsp;&nbsp;</span> <span class="OR51">&nbsp;&nbsp;&nbsp;OR51&nbsp;&nbsp;&nbsp;</span> <span class="OR52">&nbsp;&nbsp;&nbsp;OR52&nbsp;&nbsp;&nbsp;</span> <span class="OR56">&nbsp;&nbsp;&nbsp;OR56&nbsp;&nbsp;&nbsp;</span> <span class="TAAR">&nbsp;&nbsp;&nbsp;TAAR&nbsp;&nbsp;&nbsp;</span> <span class="VN1R">&nbsp;&nbsp;&nbsp;VN1R&nbsp;&nbsp;&nbsp;</span>
<span class="OR1">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR2">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR3">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR4">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR5">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR6">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR7">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR8">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR9">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR10">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR11">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR12">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR13">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR14">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR51">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR52">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="OR56">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="TAAR">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span> <span class="VN1R">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span>

<span class="rgn1">Transmembrane helix I</span>; <span class="rgn2">Transmembrane helix II</span>; <span class="rgn3">Transmembrane helix III</span>; <span class="rgn4">Transmembrane helix IV</span>; <span class="rgn5">Transmembrane helix V</span>; <span class="rgn6">Transmembrane helix VI</span>; <span class="rgn7">Transmembrane helix VII</span>; <span class="rgn8">Helix region VIII</span>.

<span class="rg50">Ballesteros-Weinstein x.50</span>; <span class="bsr">Ligand binding residue</span>; <span class="mcoord">Metal coordination residue</span>.

<span class="big5">"Big 5" receptor</span>; <i>Orphan receptor</i>; * = receptor is expressed in at least 75% of the population; ** = receptor is expressed in at least 98% of the population.
</pre>
</div>
