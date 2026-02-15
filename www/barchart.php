<?php
chdir(__DIR__);
require_once("../data/protutils.php");
require_once("../data/odorutils.php");

function lum($r, $g, $b)
{
    $rcontrib = 0.299;
    $bcontrib = 0.114;
    $gcontrib = 1.0 - $rcontrib - $bcontrib;
    return $rcontrib*$r + $gcontrib*$g + $bcontrib*$b;
}

$bkcolor = [0x15, 0x1a, 0x37];

$bkcol_tetrapod = $bkcolor;
$bkcol_tetrapod[1] += 5;

$bkcol_fishlike = $bkcolor;
$bkcol_fishlike[2] += 13;

$bkcol_taar = $bkcolor;
$bkcol_taar[0] += 10;

$bklum = lum($bkcolor[0], $bkcolor[1], $bkcolor[2]);

$odor = find_odorant(@$_REQUEST['o']);
if (!$odor)
{
    header("Location: odorants.php");
    exit;
}

$mode = 'e';            // Empirical.
if (isset($_REQUEST["m"])) $mode = $_REQUEST["m"];

$t = [];
$tprobs = [];
$e = [];
$p = [];

switch ($mode)
{
    case 'e':
    foreach ($odor['activity'] as $ref => $acv)
    {
        foreach ($acv as $rcpid => $a)
        {
            $act = @$a['adjusted_curve_top'] ?: ((@$a['type'] == 'a' && !@$a['ec50']) ? true : false);
            if (!isset($a['adjusted_curve_top']) && isset($a['type']))
            {
                $act = 0;
                if ($a['type'] == "vsa") $act = 10;
                else if ($a['type'] == "sa") $act = 8;
                else if ($a['type'] == "ma") $act = 5;
                else if ($a['type'] == "wa") $act = 2;
                else if ($a['type'] == "va") $act = 0.333;
                else if ($a['type'] == "pa")
                {
                    $act = 4;
                    $tprobs[$rcpid] = true;
                }
                if ($act) $maxt = 10;
            }

            if (!isset($t[$rcpid])) $t[$rcpid] = $act;
            else
            {
                if ($t[$rcpid] < $act) $t[$rcpid] = $act;
            }

            $ec50 = @$a['ec50'] ?: 0;
            if (!isset($e[$rcpid])) $e[$rcpid] = $ec50;
            else
            {
                if ($ec50)
                {
                    if (!$e[$rcpid] || $ec50 < $e[$rcpid]) $e[$rcpid] = $ec50;
                }
            }
        }
    }
    break;

    case 'p':
    $predictions = [];
    chdir(__DIR__);
    $dock_results = json_decode(file_get_contents("../predict/dock_results.json"), true);
    $odorname_under = str_replace(' ', '_', $odor["full_name"]);
    foreach ($dock_results as $protid => $dr)
    {
        if (isset($dr[$odorname_under]))
        {
            $p[$protid] = floatval($dr[$odorname_under]["DockScore"]);
        }
    }
    break;

    default:
    die("Unsupported mode $mode.\n");
}

$res = 3;
$xbuf = 80;
$ybuf = 80;

foreach ($prots as $protid => $prot) if (family_from_protid($protid) == "MS4A") unset($prots[$protid]);

$w = count($prots)*$res + $xbuf;
$h = 300;

$im = imagecreatetruecolor($w, $h);
imagefilledrectangle($im, 0,0, $w,$h, imagecolorallocate($im,$bkcolor[0],$bkcolor[1],$bkcolor[2]));

$last4 = "";
$xmax = count($prots);
$dxmax = $xmax * $res + $xbuf/2;
foreach (array_keys($prots) as $x => $orid)
{
    $dx  = $x * $res + $xbuf/2;

    $first4 = substr($orid, 0, 4);

    if ($first4 != $last4)
    {
        switch ($first4)
        {
            case "OR1A":
            imagefilledrectangle($im, $dx,0, $dxmax,$h, imagecolorallocate($im,$bkcol_tetrapod[0],$bkcol_tetrapod[1],$bkcol_tetrapod[2]));
            break;

            case "OR51":
            imagefilledrectangle($im, $dx,0, $dxmax,$h, imagecolorallocate($im,$bkcol_fishlike[0],$bkcol_fishlike[1],$bkcol_fishlike[2]));
            break;

            case "TAAR":
            imagefilledrectangle($im, $dx,0, $dxmax,$h, imagecolorallocate($im,$bkcol_taar[0],$bkcol_taar[1],$bkcol_taar[2]));
            break;

            default:
            ;
        }
    }

    $last4 = $first4;
}

if (!isset($maxt)) $maxt = count($t) ? ( @max($t) ?: 1 ) : 1;
$maxe = count($e) ? ( @max($e) ?: 0 ) : 0;
$mine = count($e) ? ( @min($e) ?: -6 ) : -6;
if ($maxe) $maxe += 0.5;
$mine -= 0.5;
$maxp = count($p) ? ( @max($p) ?: 1 ) : 1;

if ($maxt < 1) $maxt = 1;
if ($maxp < 1) $maxp = 1;

if ($maxe <= $mine+2) { $maxe += 1; $mine -= 1; }

if ($mine <= -3 && $maxt < 2) $maxt = 2;

$tscale = floatval($h-$ybuf) / $maxt;
$escale = floatval($h-$ybuf) / ($maxe-$mine);
$pscale = floatval($h-$ybuf) / $maxp;

$red   = imagecolorallocate($im,255,96,80);
$wine  = imagecolorallocate($im,128,32,48);
$green = imagecolorallocate($im,64,144,96);
$brown = imagecolorallocate($im,96,80,64);
$yellow= imagecolorallocate($im,192,160,96);
$blue  = imagecolorallocate($im,128,160,255);
$cyan  = imagecolorallocate($im,32,80,104);
$pink  = imagecolorallocate($im,192,176,218);
$white = imagecolorallocate($im,240,240,240);
$azure = imagecolorallocate($im,32,96,255);
$sapphire = imagecolorallocate($im,32,16,224);

function orclr($fam, $im = false)
{
    if ($im)
    {
        $rgb = orclr($fam, false);            // RECURSION!
        return imagecolorallocate($im, $rgb[0], $rgb[1], $rgb[2]);
    }

    switch ($fam)
    {
        case 1: return [0xff, 0xff, 0x99];
        case 2: return [0xff, 0xcc, 0x00];
        case 3: return [0x99, 0xff, 0x00];
        case 4: return [0x99, 0x66, 0xff];
        case 5: return [0xff, 0x66, 0x88];
        case 6: return [0x00, 0xee, 0xff];
        case 7: return [0xff, 0x33, 0x44];
        case 8: return [0xff, 0x77, 0x00];
        case 9: return [0xff, 0x55, 0x22];
        case 10: return [0x99, 0x66, 0x33];
        case 11: return [0xaa, 0xcc, 0xee];
        case 12: return [0xdd, 0xff, 0x00];
        case 13: return [0x00, 0xff, 0x99];
        case 14: return [0x33, 0x99, 0xff];
        case 51: return [0xff, 0xbb, 0x66];
        case 52: return [0xff, 0x33, 0xcc];
        case 56: return [0x22, 0xff, 0x55];
        case "TAAR": return [0x33, 0x22, 0xff];
        case "VN1R": return [0x99, 0xdd, 0xff];
        default: return [0xff, 0xff, 0xff];
    }
}

for ($i=1; $i<=14; $i++)
{
    $var = "or$i";
    $$var = orclr($i, $im);
}
for ($i=51; $i<=56; $i++)
{
    $var = "or$i";
    $$var = orclr($i, $im);
}

$taar = orclr("TAAR", $im);
$vn1r = orclr("VN1R", $im);

$base = $h-$ybuf/2;
$bsht = 8;

imageline($im, 0,$base, $w,$base, $blue);

if (count($t))
{
    imagestring($im, 3, $w-29, 0, "Rel.", $red);
    imagestring($im, 3, $w-29,15, "Top" , $red);
}
else if (count($p))
{
    imagestring($im, 3, $w-29,0, "Dock", $azure);
    imagestring($im, 3, $w-36,15, "Score" , $azure);
}

if (count($t) || count($e))
{
    // Right labels first, so that left lines take precedence.
    for ($top = 1; $top <= floor($maxt); $top += 1)
    {
        $dy = intval($base-1 - $tscale*$top);

        if (!($top & 1))
        {
            imageline($im, $xbuf/3,$dy, $w-$xbuf/3,$dy, $wine );
            imagestring($im, 3, $w-$xbuf/6,$dy-8, $top, $red);
            imagestring($im, 3, 2,$dy-8, $top, $red);
        }
    }

    // Left labels.
    /*
    imagestring($im, 3, 2,0, "log10", $green);
    imagestring($im, 3, 2,15, "EC50", $green);

    for ($ec = floor($maxe); $ec >= ceil($mine); $ec -= 1)
    {
        $dy = intval($base-1 - $escale*($maxe-$ec));

        imageline($im, $xbuf/3,$dy, $w-$xbuf/3,$dy, $cyan);
        imagestring($im, 3, 2,$dy-8, $ec, $green);
    }
        */
}

if (count($p))
{
    $step = max(0.5, floor($maxp / 7));
    for ($score = floor($maxp); $score > 0.1; $score -= $step)
    {
        $dy = intval($base-1 - $pscale*$score);

        imageline($im, $xbuf/3,$dy, $w-$xbuf/3,$dy, $sapphire);
        imagestring($im, 3, $w-$xbuf/10-$xbuf/13*strlen($score),$dy-8, $score, $azure);
    }
}

$bytree = [];
foreach ($prots as $rcpid => $pp)
{
    if (!isset($pp['btree'])) continue;
    $bytree[$pp['btree']] = $rcpid;
}

ksort($bytree, SORT_STRING);

foreach ($prots as $rcpid => $pp)
{
    if (!isset($pp['btree'])) $bytree[$rcpid] = $rcpid;
}

$texts = [];
$txtop = [];
$dybyx = [];
$cbt = count($bytree);
foreach (array_values($bytree) as $x => $orid)
{
    if (@$_REQUEST['ecexptest'])
    {
        $t[$orid] = 10;
        $e[$orid] = -(9.0 * (floatval($x)/$cbt));
        // echo $e[$orid] ."\n";
    }

    $rcpcol = $white;
    $fam = family_from_protid($orid);
    if ($fam == "MS4A") continue;
    switch (substr($orid,0,2))
    {   case 'OR':
        $bcol = 'or'.intval(preg_replace("/[^0-9]/","",substr($orid,2,2)));
        $bcol = $$bcol;
        break;

        case 'TA':
        $bcol = $taar;
        break;

        case 'VN':
        $bcol = $vn1r;
        break;

        case 'MS':
        $bcol = $ms4a;
        break;

        default:
        $bcol = $blue;
    }

    $dx  = $x * $res + $xbuf/2;
    $dyt = (isset($t[$orid]) && $t[$orid])
         ? intval($base-1 - $tscale*$t[$orid])
         : false;
    $dye = (isset($e[$orid]) && $e[$orid])
         ? intval($base-1 - $escale*($maxe-$e[$orid]))
         : false;
    $dyp = isset($p[$orid])
        ? intval($base-1 - $pscale*$p[$orid])
        : false;

    $base1 = $base2 = $base;
    if ($dyt >= $base) { $dyt += $bsht+1; $base1 += $bsht; }
    if ($dyp >= $base) { $dyp += $bsht+1; $base2 += $bsht; }

    $opacity = intval(1.1 * (100-(floatval(@$prots[$orid]['expression'] ?: 100))));
    $tred   = imagecolorallocatealpha($im,255,96,80,$opacity);
    $tgreen = imagecolorallocatealpha($im,64,144,96,$opacity);
    $tyellow= imagecolorallocatealpha($im,192,160,96,$opacity);
    $tazure = imagecolorallocatealpha($im,32,96,255,$opacity);
    if (substr($fam, 0, 2) == "OR") $orcol = orclr(intval(substr($fam, 2)));
    else $orcol = orclr($fam);

    $dy = $base;

    if (false===$dyt && false!==$dye) $dyt = $dye;

    if ($dyt)
    {
        $dy = $dyt;
        $dyscale = 1.0 / ($base-$dy);
        for ($y1 = $base; $y1 > $dyt; $y1--)
        {
            $exp = 67.1003 / pow($e[$orid] ?: 0.0001, 2);
            //if ($exp < 1) $exp = 1.0 / (1.0 - $exp);
            $opc = $e[$orid]
                ? 0.05 + 0.95 * pow(1.0-(($base-$y1) * $dyscale), $exp)
                : 0.35;
            $yc = imagecolorallocatealpha($im, $orcol[0], $orcol[1], $orcol[2], max(0, min(127, 127-127.0*$opc)));
            imagefilledrectangle($im, $dx,$y1, $e[$orid]?($dx+$res-2):$dx,$y1, $yc);
        }
    }

    /* if (false!==$dyt && false===$dye)
    {
        if (@$tprobs[$orid])
            {
                for ($y1 = $base1; $y1 > $dy=$dyt; $y1 -= 10)
                {
                    imagefilledrectangle($im, $dx,$y1, $dx+$res-2,$y1-3, $tred);
                }
            }
        else imagefilledrectangle($im, $dx,$base1, $dx+$res-2,$dy=$dyt, $tred);
    }
    if (false!==$dye && false===$dyt) imagefilledrectangle($im, $dx,$base , $dx+$res-2,$dy=$dye, $tgreen);
    if (false!==$dye && false!==$dyt)
    {
        $dy = min($dye,$dyt);
        if ($dye < $dyt)
        {
            if ($t[$orid] >= 0) imagefilledrectangle($im, $dx,$base , $dx+$res-3,$dye, $tgreen);
            imagefilledrectangle($im, $dx,$base1, $dx+$res-2,$dyt, $tyellow);
        }
        else
        {
            imagefilledrectangle($im, $dx+1,$base1, $dx+$res-2,$dyt, $tred);
            imagefilledrectangle($im, $dx,$base , $dx+$res-2,$dye, $tyellow);
        }
    } */

    if (false!==$dyp) imagefilledrectangle($im, $dx+2,$base2, $dx+$res-1,$dy=$dyp, $tazure);
    $dybyx[$x] = $dy;

    imagefilledrectangle($im, $dx,$base, $dx+$res-1,$base+$bsht, $bcol);
}

foreach (array_values($bytree) as $x => $orid)
{
    $dx  = $x * $res + $xbuf/2;
	$dy = $dybyx[$x];
    if ($dy < $h/7) $dy = $h/7;
    $peak = true;
    $xl = min($x+10, count($bytree)-1);
    for ($lx = max($x-10, 0); $lx < $xl; $lx++)
    {
    	if ($lx == $x) continue;
    	if ($dy > @$dybyx[$lx]) $peak = false;
    }
    if ($peak && $dy < $h/1.4)
    {
        $texts[] = [$dx, $dy-5, @$tprobs[$orid] ? "($orid)" : $orid];
        $txtop[] = intval(0.81 * (100-(floatval(@$prots[$orid]['expression'] ?: 100))));
    }
}

if (!file_exists($fontfile = "assets/Montserrat.ttf"))
{
    $wofffile = str_replace(".ttf", ".woff2", $fontfile);
    $c = file_get_contents("https://fonts.gstatic.com/s/montserrat/v25/JTUHjIg1_i6t8kCHKm4532VJOt5-QNFgpCtr6Hw5aXo.woff2");
    $fp = fopen($wofffile, "wb");
    if (!$fp) die("Failed to open $wofffile for writing!");
    fwrite($fp, $c);
    fclose($fp);

    exec("woff2_decompress $wofffile");
}


$x = 53*$res;
$x1 = $x + 10*strlen($odor['full_name']);

foreach ($texts as $k => $txt)
{
    $tpink = imagecolorallocatealpha($im,192,176,218,$txtop[$k]);
    imagettftext($im, 9, 35, $txt[0], $txt[1], $tpink, $fontfile, $txt[2]);

    if ($txt[0] >= $x && $txt[0] <= $x1) $x = $txt[0] + 50;
}

if ((!count($t) || !max($t)) && (!count($e) || !min($e)) && !count($p))
{   
    imagettftext($im, 28, 0, $w*0.42, $h*0.44, $pink, $fontfile, "(no data)");
}

// Odorant Name
imagestring($im, 5, $x, 2, $odor['full_name'], $blue);

header('Content-Type: image/png');
header('Content-Disposition: inline; filename="barchart.png"');
imagepng($im);
