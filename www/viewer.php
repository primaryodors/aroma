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
    $path = "../out/$fam/$protid/$protid.$odor.$mode.model$n.pdb";
    if (!file_exists($path)) $path = "../out/$fam/$protid/$protid.$odor.model$n.pdb";
    if (!file_exists($path)) die("Something went wrong.");
    $pdb = file_get_contents($path);

    $dock = "../out/$fam/$protid/$protid.$odor.$mode.dock";

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
    require("../data/3dtools.php");
    require("../data/protutils.php");
    require("../data/odorutils.php");
    $protid = $_REQUEST["prot"];
    $fam = family_from_protid($protid);
    $sub = subfamily_from_protid($protid);
    $odor = $_REQUEST["odor"];
    $mode = $_REQUEST["mode"];      // active or inactive.
    $n = @$_REQUEST["n"] ?: 1;

    // SVG section parameters and variables
    $wid = 602; $hei = 420; $scale = floatval($wid) / 20;
    $resradius = 25;
    $capos = [];
    $ligapos = [];
    $ligcen = [0,0,0];
    $ligbonds = [];
    $farthest1 = $farthest2 = false;
    $farthestr = 0;
    $ligrot1 = [0,0,0,0];
    $ligrot2 = [0,0,0,0];
    $ligrot3 = [0,0,0,0];
    $bsrintera = [];

    chdir(__DIR__);
    $dock = "../out/$fam/$protid/$protid~$odor.$mode.dock";
    if (!file_exists($dock)) die("Something went wrong.");
    $txt = file_get_contents($dock);

    $cavfn = "../pdbs/$fam/$protid.$mode.cvty";
    if (!file_exists($cavfn)) $cavfn = "../cpllocal/$fam/$sub/$protid.$mode.cvty";
    if (file_exists($cavfn))
    {
        $lines = explode("\n", $txt);
        $poseno = 0;
        foreach ($lines as $i => $ln)
        {
            if (substr($ln, 0, 6) == "Pose: ")
            {
                $poseno = intval(substr($ln, 6));
            }
            else if ($poseno == 1 && substr($ln, 0, 4) == "ATOM")
            {
                $aname = trim(substr($ln, 12, 4));
                if ($aname == "CA")
                {
                    $resno = intval(substr($ln, 22, 4));
                    $x = floatval(substr($ln, 30, 8));
                    $y = floatval(substr($ln, 38, 8));
                    $z = floatval(substr($ln, 46, 8));
                    $capos[$resno] = [$x,$y,$z];
                }
            }
            else if ($poseno == 1 && substr($ln, 0, 6) == "HETATM")
            {
                $aname = trim(substr($ln, 12, 4));
                $x = floatval(substr($ln, 30, 8));
                $y = floatval(substr($ln, 38, 8));
                $z = floatval(substr($ln, 46, 8));
                $ligapos[$aname] = [$x,$y,$z];
            }
            else if ($poseno == 1 && preg_match("/^[A-Z][a-z]{2}[0-9]+:[A-Z0-9]+~[^:]+:[A-Za-z0-9]+:\\s+[0-9.+-]+/", $ln))
            {
                list($resno, $idgaf, $aname, $binding) = explode(':', $ln);
                if (false===strpos($binding, ' ')) $binding .= " none";
                list($energy, $intera_type) = explode(' ', trim($binding));
                $resno = intval(preg_replace("/[^0-9]/", "", $resno));
                $energy = floatval($energy);

                if ($energy < -0.5)
                {
                    $bw = bw_from_resno($protid, $resno);
                    $bsrintera[$bw][$aname] = [ $energy, $intera_type ];
                }
            }
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

    if ($nla = count($ligapos))
    {
        foreach ($ligapos as $aname => $xyz)
        {
            $ligcen[0] += $xyz[0];
            $ligcen[1] += $xyz[1];
            $ligcen[2] += $xyz[2];

            if (substr($aname, 0, 1) == "H") continue;
            foreach ($ligapos as $bname => $abc)
            {
                if ($bname == $aname) continue;
                if (substr($bname, 0, 1) == "H") continue;
                $r = sqrt( pow($xyz[0] - $abc[0], 2)
                         + pow($xyz[1] - $abc[1], 2)
                         + pow($xyz[2] - $abc[2], 2)
                         );
                if ($r < 2.2)
                    $ligbonds[$aname][$bname] = true;
                if ($r > $farthestr)
                {
                    $farthest1 = $aname;
                    $farthest2 = $bname;
                    $farthestr = $r;
                }
            }
        }

        foreach ($ligcen as $k => $v) $ligcen[$k] /= $nla;
    }

    if ($farthestr) $ligrot1 = align_points_3d($ligapos[$farthest2], [10000,$ligcen[1],$ligcen[2]], $ligcen);

    $ligaxy = [];
    $fiftyseventh = pi()/180;
    $bestspread = 0.0;
    $besttheta = 0.0;
    for ($theta=0.0; $theta<pi()*2; $theta+=5.0*$fiftyseventh)
    {
        $ligrot2 = [ 10000.0, 0.0, 0.0, $theta ];
        foreach ($ligapos as $aname => $xyz)
        {
            // 3D rotation
            $xyz = rotate3D($xyz, $ligcen, $ligrot1, $ligrot1[3]);
            $xyz = rotate3D($xyz, $ligcen, $ligrot2, $ligrot2[3]);
            list($x, $y, $z) = $xyz;

            $cx = intval($wid/2 + ($x - $ligcen[0])*$scale);
            $cy = intval($hei/2 + ($z - $ligcen[2])*$scale);
            $ligaxy[$aname] = [$cx, $cy];
        }

        $lspread = 0;
        foreach ($ligaxy as $aname => $cxy)
        {
            foreach ($ligaxy as $bname => $dxy)
            {
                if ($bname == $aname) continue;
                $r = get_2d_distance($cxy, $dxy);
                $lspread += $r;
            }
        }

        if ($lspread > $bestspread)
        {
            $bestspread = $lspread;
            $besttheta = $theta;
        }
    }

    $ligrot2 = [ 10000.0, 0.0, 0.0, $besttheta ];
    foreach ($ligapos as $aname => $xyz)
    {
        // 3D rotation
        $xyz = rotate3D($xyz, $ligcen, $ligrot1, $ligrot1[3]);
        $xyz = rotate3D($xyz, $ligcen, $ligrot2, $ligrot2[3]);
        list($x, $y, $z) = $xyz;

        $cx = intval($wid/2 + ($x - $ligcen[0])*$scale);
        $cy = intval($hei/2 + ($z - $ligcen[2])*$scale);
        $ligaxy[$aname] = [$cx, $cy];
    }

    $besttheta = 0.0;
    $bestspread = 1e9;
    for ($theta=0.0; $theta<pi()*2; $theta+=5.0*$fiftyseventh)
    {
        $ligrot3 = [0, 0, 10000, $theta];
        $ymin = 1e9;
        $ymax = -1e9;
        foreach ($ligaxy as $aname => $cxy)
        {
            // 2D rotation trick
            $cxy[2] = 0;
            $cxy = rotate3D($cxy, $ligcen, $ligrot3, $ligrot3[3]);
            list($cx, $cy, $cz) = $cxy;

            if ($cy < $ymin) $ymin = $cy;
            if ($cy > $ymax) $ymax = $cy;
        }

        $yspread = $ymax - $ymin;
        if ($yspread < $bestspread)
        {
            $bestspread = $yspread;
            $besttheta = $theta;
        }
    }
    $ligrot3 = [0, 0, 10000, $besttheta];
    foreach ($ligapos as $aname => $xyz)
    {
        // 3D rotation
        $xyz = rotate3D($xyz, $ligcen, $ligrot1, $ligrot1[3]);
        $xyz = rotate3D($xyz, $ligcen, $ligrot2, $ligrot2[3]);
        list($x, $y, $z) = $xyz;

        $cx = intval($wid/2 + ($x - $ligcen[0])*$scale);
        $cy = intval($hei/2 + ($z - $ligcen[2])*$scale);
        $ligaxy[$aname] = rotate3D([$cx, $cy, 0], [$wid/2, $hei/2, 0], $ligrot3, $ligrot3[3]);
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
    if (isset($gprots["hGNA$mode"])) echo "<h1>$protid:GNA$mode ~ $odor</h1>";
    else echo "<h1>$protid ~ $odor</h1>";
    $bsr4sim = array_values(similar_receptors($protid))[0][0];
    // print_r($bsr4sim);

    $dockfname = "../out/$fam/$protid/$protid~$odor.$mode.dock";
    $lbsr = [];
    $lbstr = [];
    $d = file_get_contents($dockfname);
    $lines = explode("\n", $d);
    foreach ($lines as $ln) 
    {
        if (trim($ln) == "# PDB Data") break;
        if (trim($ln) == "Pose: 2") break;
        if (false!==strpos($ln, "~(ligand)") || false!==strpos($ln, "$odor:"))
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
        try
        {
            $lbsrn[$bw] = resno_from_bw($protid, $bw);
        }
        catch (Exception $ex)
        {
            continue;
        }
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
        foreach ($sim as $id => list($lb, $sim_percent))
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
            $sim_percent = intval($sim_percent*100);
            echo "$sim_percent%";
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
                    var svgdat = "<svg id=\"function random() { [native code] }\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"<?php echo $wid; ?>px\" height=\"<?php echo $hei; ?>px\" viewBox=\"0 0 <?php echo $wid; ?> <?php echo $hei; ?>\">"; // svg_from_smiles("<?php echo $o["smiles"]; ?>", <?php echo $wid; ?>, <?php echo $hei; ?>);
                    svgdat = svgdat.replace("</svg>", "");
                    <?php

                    // RESIDUE RAW X,Y COORDINATES FOR SVG
                    $rescxy = [];
                    foreach ($lb as $bw => $aa)
                    {
                        $i = intval(preg_replace("/[^0-9]/", "", $lbsr[$bw])) - 1;
                        if ($i < 0) continue;

                        $resno = resno_from_bw($protid, $bw);
                        if (@$caloc = $capos[$resno])
                        {
                            // 3D rotation
                            $lxyz = rotate3D($caloc, $ligcen, $ligrot1, $ligrot1[3]);
                            list($x, $y, $z) = rotate3D($lxyz, $ligcen, $ligrot2, $ligrot2[3]);

                            $cx = intval($wid/2 + ($x - $ligcen[0])*$scale);
                            $cy = intval($hei/2 + ($z - $ligcen[2])*$scale);
                            $rescxy[$bw] = rotate3D([$cx, $cy, 0], [$wid/2, $hei/2, 0], $ligrot3, $ligrot3[3]);
                        }
                    }

                    // MOVE RESIDUES TO TIDY LOCATIONS
                    for ($respositer=0; $respositer<503; $respositer++)
                    {
                        foreach ($rescxy as $bw => $cxy)
                        {
                            $nessamo = false;
                            $nessamd = 1e9;
                            $nessamxy = [0,0];

                            foreach ($ligaxy as $aname => $lxy)
                            {
                                $r = get_2d_distance($cxy, $lxy);
                                if ($r < $nessamd)
                                {
                                    $nessamd = $r;
                                    $nessamo = $aname;
                                    $nessamxy = $lxy;
                                }
                            }

                            if ($nessamo)
                            {
                                $r = get_2d_distance($nessamxy, $cxy);
                                $adjustment = [ $nessamxy[0] - $cxy[0], $nessamxy[1] - $cxy[1] ];
                                $d = $nessamd - 90;
                                $adjustment[0] *= $d / $r / 3;
                                $adjustment[1] *= $d / $r / 3;

                                $rescxy[$bw] = [ $cxy[0]+$adjustment[0], $cxy[1]+$adjustment[1] ];
                            }

                            foreach ($rescxy as $ne => $vcu)
                            {
                                if ($ne == $bw) continue;
                                $r = get_2d_distance($cxy, $vcu);
                                if ($r < 65)
                                {
                                    $adjustment = [ ($vcu[0]-$cxy[0])*10.0/$r, ($vcu[1]-$cxy[1])*10.0/$r ];
                                    $rescxy[$bw] = [ $cxy[0]-$adjustment[0], $cxy[1]-$adjustment[1] ];
                                    $rescxy[$ne] = [ $vcu[0]+$adjustment[0], $vcu[1]+$adjustment[1] ];
                                }
                            }
                        }
                    }

                    // INTERACTION LINES
                    foreach ($bsrintera as $bw => $lresinters)
                    {
                        foreach ($lresinters as $aname => list($aname, $intera_type))
                        {
                            if (isset($ligaxy[$aname]) && isset($rescxy[$bw]))
                            {
                                $dasharray = "stroke-dasharray=\\\"2\\\"";
                                switch ($intera_type)
                                {
                                    case "hbond": case "ionic":
                                    $couleur = "#393";
                                    break;

                                    case "polarpi": case "mcoord":
                                    $couleur = "#696";
                                    break;

                                    case "pi":
                                    $couleur = "#c9c";
                                    break;

                                    case "covalent":
                                    $couleur = "#ccc";
                                    $dasharray = "";
                                    break;

                                    default:
                                    $couleur = "#666";
                                }

                                list($x1, $y1) = $ligaxy[$aname];
                                list($x2, $y2) = $rescxy[$bw];
                                $r = get_2d_distance($ligaxy[$aname], $rescxy[$bw]);
                                $coeff = ($r-$resradius)/$r;
                                $x2 = $x1 + $coeff*($x2-$x1);
                                $y2 = $y1 + $coeff*($y2-$y1);

                                echo "svgdat += \"<line x1=\\\"$x1\\\" y1=\\\"$y1\\\" x2=\\\"$x2\\\" y2=\\\"$y2\\\" stroke=\\\"$couleur\\\" $dasharray />\\n\";\n";
                            }
                        }
                    }

                    // LIGAND COVALENT BONDS
                    foreach ($ligbonds as $aname1 => $b2)
                    {
                        if (!isset($ligaxy[$aname1])) continue;
                        list($x1,$y1) = $ligaxy[$aname1];
                        foreach ($b2 as $aname2 => $v)
                        {
                            if ($v && isset($ligaxy[$aname2]))
                            {
                                list($x2,$y2) = $ligaxy[$aname2];
                                echo "svgdat += \"<line x1=\\\"$x1\\\" y1=\\\"$y1\\\" x2=\\\"$x2\\\" y2=\\\"$y2\\\" stroke=\\\"#ccc\\\" />\\n\";\n";
                            }
                        }
                    }

                    // LIGAND ATOMS
                    foreach ($ligapos as $aname => $xyz)
                    {
                        $elem = ucwords(strtolower(preg_replace("/[0-9]/", "", $aname)));
                        if ($elem == "H") continue;
                        switch ($elem)
                        {
                            case "C":
                                $couleur = "#ccc";
                                break;
                            case "N":
                                $couleur = "#06f";
                                break;
                            case "O":
                                $couleur = "#f00";
                                break;
                            case "S":
                                $couleur = "#fc0";
                                break;
                            case "Cl":
                                $couleur = "#0c0";
                                break;
                            case "Br":
                                $couleur = "#c60";
                                break;
                            case "I":
                                $couleur = "#609";
                                break;
                            default:
                                $couleur = "#f6f";
                                break;
                        }
                        list($cx,$cy) = $ligaxy[$aname];
                        echo "svgdat += \"<circle cx=\\\"$cx\\\" cy=\\\"$cy\\\" r=\\\"4\\\" fill=\\\"$couleur\\\"></circle>\\n\";\n";
                    }

                    // RESIDUE CIRCLES
                    foreach ($lb as $bw => $aa)
                    {
                        if (isset($rescxy[$bw]))
                        {
                            $resno = resno_from_bw($protid, $bw);
                            list($cx,$cy) = $rescxy[$bw];
                            $tx  = $cx - intval(floatval(strlen($bw)) * 3.333);
                            $txr = $cx - intval(floatval(strlen($resno)+1) * 3.333);
                            $ty  = $cy - 2;
                            $tyr = $cy + 9;

                            switch ($aa)
                            {
                                case "A": case "I": case "L": case "V": case "P":
                                    $couleur = "#ccc";
                                    break;
                                case "M":
                                    $couleur = "#cc6";
                                    break;
                                case "G":
                                    $couleur = "#666";
                                    break;
                                case "C":
                                    $couleur = "#fc0";
                                    break;
                                case "F": case "W":
                                    $couleur = "#c9c";
                                    break;
                                case "Y":
                                    $couleur = "#f6f";
                                    break;
                                case "S": case "T": case "N": case "Q":
                                    $couleur = "#0ff";
                                    break;
                                case "D": case "E":
                                    $couleur = "#f66";
                                    break;
                                case "K": case "R":
                                    $couleur = "#06f";
                                    break;
                                case "H":
                                    $couleur = "#76f";
                                    break;
                                default:
                                    $couleur = "#f0f";
                            }

                            echo "svgdat += \"<circle cx=\\\"$cx\\\" cy=\\\"$cy\\\" r=\\\"$resradius\\\" stroke=\\\"$couleur\\\" stroke-width=\\\"1\\\" fill-opacity=\\\"0\\\"></circle>\\n\";\n";
                            echo "svgdat += \"<text x=\\\"$tx\\\" y=\\\"$ty\\\" fill=\\\"$couleur\\\" font-size=\\\"11px\\\">$bw</text>\\n\";\n";
                            echo "svgdat += \"<text x=\\\"$txr\\\" y=\\\"$tyr\\\" fill=\\\"$couleur\\\" font-size=\\\"11px\\\">$aa$resno</text>\\n\";\n";
                        }
                    }
                    ?>
                    svgdat += "</svg>";
                    $("#skel<?php echo $id; ?>")[0].innerHTML = svgdat;
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
                    if (0) foreach ($lb as $bw => $aa)
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
                    $bwx = str_replace('.', 'x', $bw);
                    $display = isset($lbsr[$bw]) ? "" : "display: none;";
                    echo "<th style=\"$display\" class=\"show$bwx\">$bw</th>\n";
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
            <td colspan="2" class="bigger">tighter</td>
            <td colspan="2" class="smaller">roomier</td>
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
