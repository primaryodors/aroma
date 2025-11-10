<?php
chdir(__DIR__);

$c = file_get_contents("../viewer.htm");

if (@$_REQUEST['view'] == "pred")
{
    require("../data/protutils.php");
    $protid = $_REQUEST["prot"];
    $fam = family_from_protid($protid);
    $odor = $_REQUEST["odor"];
    $mode = $_REQUEST["mode"];      // active or inactive.
    $n = @$_REQUEST["n"] ?: 1;

    chdir(__DIR__);
    $path = "../output/$fam/$protid/$protid.$odor.$mode.model$n.pdb";
    if (!file_exists($path)) $path = "../output/$fam/$protid/$protid.$odor.model$n.pdb";
    if (!file_exists($path)) die("Something went wrong.");
    $pdb = file_get_contents($path);

    $dock = "../output/$fam/$protid/$protid.$odor.$mode.dock";

    $ligbs = "	var lligbs = [";
    $pdblines = explode("\n", $pdb);
    $atom_lines = [];
    $hetatm_xyz = [];
    foreach ($pdblines as $line)
    {
        if (substr($line, 0, 7) == "ATOM   ") $atom_lines[] = $line;
        if (substr($line, 0, 7) == "HETATM ")
        {
            $x = floatval(substr($line, 31, 7));
            $y = floatval(substr($line, 39, 7));
            $z = floatval(substr($line, 47, 7));
            $hetatm_xyz[] = [$x, $y, $z];
        }
    }

    $ligand_residues = [];
    foreach ($atom_lines as $line)
    {
        $resno = intval(substr($line, 23, 3));
        if (isset($ligand_residues[$resno])) continue;
        $x = floatval(substr($line, 31, 7));
        $y = floatval(substr($line, 39, 7));
        $z = floatval(substr($line, 47, 7));

        foreach ($hetatm_xyz as $xyz)
        {
            $r = sqrt(pow($x-$xyz[0], 2) + pow($y-$xyz[1], 2) + pow($z-$xyz[2], 2));
            if ($r <= 5.0)
            {
                if (count($ligand_residues)) $ligbs .= ", ";
                $ligbs .= "$resno";
                $ligand_residues[$resno] = $resno;
                break;
            }
        }
    }
    $ligbs .= "];\n";

    $c = str_replace("	var lligbs = get_ligbs_from_orid();\n", $ligbs, $c);
    $c = str_replace("var literal_pdb = false;\n", "var literal_pdb = `$pdb`;\n", $c);
    $c = str_replace("var literal_fname = \"\";\n", "var literal_fname = \"$protid.$odor.$mode.model$n.pdb\";\n", $c);

    $d = explode("\n", file_get_contents($dock));
    $dockdisp = [];
    $ddidx = 0;
    $dockdisp[$ddidx] = "";
    foreach ($d as $lineno => $line)
    {
        $line2 = trim($d[$lineno+2]);
        if ($line2 == "PDBDAT:") break;

        $dockdisp[$ddidx] .= "$line\n";
    }

    $c .= <<<dockdata

<script>
$('#dockfloat span')[0].innerText = `{$dockdisp[0]}`;
</script>
dockdata;

}

if (@$_REQUEST['view'] == "dock")
{
    require("../data/protutils.php");
    $protid = $_REQUEST["prot"];
    $fam = family_from_protid($protid);
    $odor = $_REQUEST["odor"];
    $mode = $_REQUEST["mode"];      // active or inactive.
    $n = @$_REQUEST["n"] ?: 1;

    chdir(__DIR__);
    $dock = "../output/$fam/$protid/$protid~$odor.$mode.dock";
    if (!file_exists($dock)) die("Something went wrong.");
    $txt = file_get_contents($dock);

    $cavfn = "../pdbs/$fam/$protid.".($mode=='active'?$mode:"upright").".cvty";
    if (file_exists($cavfn))
    {
        $lines = explode("\n", $txt);
        foreach ($lines as $i => $ln)
        {
            if (substr($ln, 0, 10) != "REMARK 800") continue;
            $next = @$lines[$i+1];
            if (substr($next, 0, 10) != "REMARK 800")
            {
                $ln .= "\nREMARK 821";
                $cavs = explode("\n", file_get_contents($cavfn));
                foreach ($cavs as $cav)
                {
                    if (!trim($cav)) continue;
                    $ln .= "\nREMARK 821 $cav";
                }
                $ln .= "\nREMARK 821";
                $lines[$i] = $ln;
                break;
            }
        }
    
        $txt = implode("\n", $lines);
    }

    // $c = str_replace("	var lligbs = get_ligbs_from_orid();\n", $ligbs, $c);
    $c = str_replace("var literal_pdb = false;\n", "var literal_pdb = `$txt`;\n", $c);
    $c = str_replace("var literal_fname = \"\";\n", "var literal_fname = \"$protid~$odor.$mode.dock\";\n", $c);

    $c .= <<<dockdata

<script>
$('#dockfloat span')[0].innerText = `{$dockdisp[0]}`;
</script>
dockdata;

    echo "<link rel=\"stylesheet\" href=\"assets/style.css\">\n";
    echo "<div style=\"display: flex;\">";
    echo "<div style=\"display: block; width: 50%; height: 100vh; overflow: auto;\">";

    echo "<h1>$protid ~ $odor</h1>";

    $dockfname = "../output/$fam/$protid/$protid~$odor.$mode.dock";
    // if (!file_exists($dockfname)) die("oops");
    $lbsr = [];
    $lbstr = [];
    $d = file_get_contents($dockfname);
    $lines = explode("\n", $d);
    foreach ($lines as $ln) 
    {
        if (trim($ln) == "# PDB Data") break;
        if (trim($ln) == "Pose: 2") break;
        if (false!==strpos($ln, "~(ligand)"))
        {
            list($lcntct, $strength) = explode(": ", $ln, 2);
            $strength = floatval($strength);
            if ($strength <= -0.1)
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
                if (!isset($lbstr[$bw]) || $strength < $lbstr[$bw])
                {
                    $lbsr[$bw] = $liga;
                    $lbstr[$bw] = $strength;
                }
            }
        }
    }
    $sim = similar_receptors($protid); // , array_keys($lbsr));
    $lbsrn = [];
    $frist = true;
    foreach (array_keys($sim[$protid]) as $bw)
    {
        if (!$frist) echo " | ";
        $lbsrn[$bw] = resno_from_bw($protid, $bw);
        $bwx = str_replace('.', 'x', $bw);
        echo "<a href=\"#\" onclick=\"$('.show$bwx').toggle();\">{$lbsrn[$bw]}</a>";
        $frist = false;
    }
    ?>
    <table class="simr">
        <?php
        $frist = true;
        foreach ($sim as $id => $lb)
        {
            if ($frist)
            {
                echo "<tr>\n";
                echo "<th>&nbsp;</th>";
                foreach (array_keys($lb) as $bw)
                {
                    $display = isset($lbsr[$bw]) ? "" : "display: none;";
                    $resno = $lbsrn[$bw];
                    $bwx = str_replace('.', 'x', $bw);
                    echo "<th style=\"$display\" class=\"show$bwx\">$resno</th>\n";
                }
                echo "</tr>\n";
                echo "<tr>\n";
                echo "<th>&nbsp;</th>";
                foreach (array_keys($lb) as $bw)
                {
                    $display = isset($lbsr[$bw]) ? "" : "display: none;";
                    $bwx = str_replace('.', 'x', $bw);
                    echo "<th style=\"$display\" class=\"show$bwx\">$bw</th>\n";
                }
                echo "</tr>\n";
            }
            echo "<tr>\n";
            echo "<th><a href=\"receptor.php?r=$id\">$id</a></th>";
            foreach ($lb as $bw => $aa)
            {
                $display = isset($lbsr[$bw]) ? "" : "display: none;";
                $bwx = str_replace('.', 'x', $bw);
                echo "<td class=\"aacolor$aa show$bwx\" style=\"$display\">$aa</td>\n";
            }
            echo "</tr>\n";
            if ($frist)
            {
                echo "<tr>\n";
                echo "<th>&nbsp;</th>";
                foreach (array_keys($lb) as $bw)
                {
                    $display = isset($lbsr[$bw]) ? "" : "display: none;";
                    $la = $lbsr[$bw];
                    $bwx = str_replace('.', 'x', $bw);
                    echo "<th style=\"$display\" class=\"show$bwx\">$la</th>\n";
                }
                echo "</tr>\n";
            }
            $frist = false;
        }
        ?>
    </table>
    <?php

    echo "</div>";
    echo "<div style=\"display: block; width: 50%; height: 100vh;\">";
    $c = str_replace("<div id=\"viewport\" style=\"width:1800px;", "<div id=\"viewport\" style=\"width:60%;", $c);
    echo $c;
    echo "</div>";
    echo "</div>";
    exit;
}

echo $c;
