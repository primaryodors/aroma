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
    require("../data/odorutils.php");
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
    $c = preg_replace("/<title>[A-Za-z0-9 -]+</", "<title>$protid~$odor<", $c);

    $c .= <<<dockdata

<script>

$('#dockfloat span')[0].innerText = `{$dockdisp[0]}`;
</script>
dockdata;

    echo "<link rel=\"stylesheet\" href=\"assets/style.css\">\n";
    echo "<div style=\"display: flex;\">";
    echo "<div style=\"display: block; width: 50%; height: 100vh; padding-left: 15px; overflow: auto;\">";

    ?>
<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.6.0/jquery.min.js"></script>
<script src="https://www.lactame.com/lib/openchemlib/5.2.0/openchemlib-minimal.js"></script>
<script>

var ax = [], ay = [];
function svg_from_smiles(smiles, w, h)
{
    var molecule=OCL.Molecule.fromSmiles(smiles);
    result = molecule.toSVG(w, h, Math.random.toString(36), {fontWeight: 900})
        .replace(/rgb\(0,0,0\)/g,"rgb(255,255,255)")
        .replace(/fill=\"rgb\(160,0,0\)\">.*<\/text/g, '></text')
        .replace(/rgb\(160,0,0\)/g,"rgb(170,187,204)")
        ;
    var n = molecule.getAllAtoms();
    var i;
    for (i=0; i<n; i++)
    {
        ax[i] = molecule.getAtomX(i);
        ay[i] = molecule.getAtomY(i);
        // console.log(molecule.getAtomLabel(i) + ": X="+ax[i] + ", Y="+ay[i]);
    }
    return result;
}
</script>
<?php
    echo "<h1>$protid ~ $odor</h1>";
    $bsr4sim = array_values(similar_receptors($protid))[0][0];
    // print_r($bsr4sim);

    $dockfname = "../output/$fam/$protid/$protid~$odor.$mode.dock";
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
                if (!isset($bsr4sim[$bw])) $bsr4sim[$bw] = $resa;
            }
        }
    }

    ksort($bsr4sim);
    // print_r($bsr4sim);

    $sim = similar_receptors($protid, array_keys($bsr4sim));
    $lbsrn = [];
    echo "<p>Toggle:";
    // print_r($sim);
    foreach (array_keys($sim[$protid][0]) as $bw)
    {
        echo " &nbsp; ";
        $lbsrn[$bw] = resno_from_bw($protid, $bw);
        $aa = letter_at_bw($protid, $bw);
        $bwx = str_replace('.', 'x', $bw);
        echo "<a href=\"#\" onclick=\"$('.show$bwx').toggle();\">$aa{$lbsrn[$bw]}</a>";
    }
    echo "</p>";
    ?>
    <table class="simr">
        <?php
        $frist = true;
        $lataa = [];
        $o = find_odorant($odor);
        foreach ($sim as $id => list($lb, $simpcnt))
        {
            if ($fam != family_from_protid($id)) continue;
            if ($frist)
            {
                echo "<tr>\n";
                echo "<th>&nbsp;</th>";
                echo "<th>&nbsp;</th>";
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
                echo "<th>Prot.</th>";
                echo "<th>Sim.</th>";
                echo "<th>Expr.</th>";
                foreach (array_keys($lb) as $bw)
                {
                    $display = isset($lbsr[$bw]) ? "" : "display: none;";
                    $bwx = str_replace('.', 'x', $bw);
                    echo "<th style=\"$display\" class=\"show$bwx\">$bw</th>\n";
                }
                echo "</tr>\n";
            }
            echo "<tr onclicqk=\"$('.skeltr').hide(); $('#skeltr$id').show();\">\n";
            echo "<td style=\"text-align: left;\">";
            echo "<b><a href=\"receptor.php?r=$id\">$id</a></b>";
            echo "</td>";
            echo "<td style=\"text-align: left;\">";
            $simpcnt = intval($simpcnt*100);
            echo "$simpcnt%";
            echo "</td>";
            echo "<td style=\"text-align: left;\">";
            if (isset($prots[$id]['expression'])) echo " {$prots[$id]['expression']}%";
            echo "</td>";
            foreach ($lb as $bw => $aa)
            {
                $display = isset($lbsr[$bw]) ? "" : "display: none;";
                $bwx = str_replace('.', 'x', $bw);
                if ($frist) $lataa[$bw] = $aa;
                $relstyle = "";
                if ($lataa[$bw] != $aa)
                {
                    $cmp = amino_cmp_size($aa, $lataa[$bw]);
                    if ($cmp > 0) $relstyle = "bigger";
                    else if ($cmp < 0) $relstyle = "smaller";
                    else
                    {
                        $cmp = amino_cmp_pi($aa, $lataa[$bw]);
                        if ($cmp > 0) $relstyle = "ener";
                        else if ($cmp < 0) $relstyle = "phater";
                    }

                    $cmp = amino_cmp_hydro($aa, $lataa[$bw]);
                    if ($cmp > 0) $relstyle = "wetter";
                    else if ($cmp < 0) $relstyle = "drier";

                    $cmp = amino_cmp_charge($aa, $lataa[$bw]);
                    if ($cmp > 0) $relstyle = "basic";
                    else if ($cmp < 0) $relstyle = "acidic";
                }
                echo "<td class=\"aacolor$aa show$bwx $relstyle\" style=\"$display\">$aa</td>\n";
            }
            echo "</tr>\n";
            echo "<tr class=\"skeltr\" id=\"skeltr$id\"";
            if (!$frist) echo " style=\"display: none;\"";
            echo ">\n";
            if ($frist)
            {
                ?>
                <td colspan="16">
                <div id="skel<?php echo $id; ?>"></div>
                </td>
                </tr>
                <script>
                window.setTimeout(function()
                {
                    $("#skel<?php echo $id; ?>")[0].innerHTML = svg_from_smiles("<?php echo $o["smiles"]; ?>", 300, 300);
                    var cx = 0, cy = 0, x0, x1, y0, y1;
                    var i, n = ax.length;
                    for (i=0; i<n; i++)
                    {
                        if (!i || ax[i] < x0) x0 = ax[i];
                        if (!i || ax[i] > x1) x1 = ax[i];
                        if (!i || ay[i] < y0) y0 = ay[i];
                        if (!i || ay[i] > y1) y1 = ay[i];
                    }
                    cx = (x0 + x1) / 2;
                    cy = (y0 + y1) / 2;
                    var rect = $("#skel<?php echo $id; ?> svg")[0].getClientRects()[0];
                    console.log(rect);
                    <?php
                    $aayoff = [];
                    foreach ($lb as $bw => $aa)
                    {
                        $i = intval(preg_replace("/[^0-9]/", "", $lbsr[$bw])) - 1;
                        $ay = @$aayoff[$i] ?: 0;
                        if ($i < 0) continue;
                        echo "                var x = parseInt((ax[$i]-cx)*40-8), y = parseInt((ay[$i]-cy)*50+13*$ay), ih = '$aa<sup>$bw</sup>', cls = 'aacolor$aa';\n"; ?>
                        var aa = document.createElement("span");
                        aa.innerHTML = ih;
                        aa.className = cls;
                        aa.style.position = 'absolute';
                        if (x < 0) x -= 15;
                        if (y < 20) y -= 20;
                        aa.style.top = parseInt(rect.top + rect.height/2 + y).toString() + "px";
                        aa.style.left = parseInt(rect.left + rect.width/2 + x).toString() + "px";
                        $("#skel<?php echo $id; ?>")[0].appendChild(aa);
                        <?php
                        $aayoff[$i] = $ay+1;
                    }
                    ?>
                }, 250);
                </script>
                <?php
            }
            if ($frist)
            {
                echo "<tr>\n";
                echo "<th colspan=\"3\">&nbsp;</th>";
                foreach (array_keys($lb) as $bw)
                {
                    $display = isset($lbsr[$bw]) ? "" : "display: none;";
                    echo "<th style=\"$display\" class=\"show$bwx\">$bw</th>\n";
                    /* $la = $lbsr[$bw];
                    $bwx = str_replace('.', 'x', $bw);
                    echo "<th style=\"$display\" class=\"show$bwx\">$la</th>\n"; */
                }
                echo "</tr>\n";
            }
            $frist = false;
        }
        ?>
        <tr>
            <td colspan="16">&nbsp;</td>
        </tr>
        <tr>
            <th>Legend:</th>
            <td colspan="2" class="bigger">bigger</td>
            <td colspan="2" class="smaller">smaller</td>
            <td colspan="2" class="wetter">wetter</td>
            <td colspan="2" class="drier">drier</td>
            <td colspan="2" class="acidic">more acid</td>
            <td colspan="2" class="basic">more base</td>
            <td colspan="2" class="ener">more pi</td>
            <td colspan="2" class="phater">less pi</td>
        </tr>
    </table>
    <?php

    echo "</div>";
    echo "<div id=\"wrapport\" style=\"background-color: #000; display: block; width: 50%; height: 100vh;\">";
    $c = str_replace("<div id=\"viewport\" style=\"width:1800px;", "<div id=\"viewport\" style=\"width:60%;", $c);
    echo $c;
    echo "</div>";
    echo "</div>";
    exit;
}

echo $c;
