<?php
header("Access-Control-Allow-Origin: *");

header("Cache-Control: no-cache, no-store, must-revalidate");
header("Pragma: no-cache");
header("Expires: 0");

chdir(__DIR__);
require_once("../data/odorutils.php");

if (@$_REQUEST['m'] == "rand") $odor = $odors[array_keys($odors)[rand(0,count($odors)-1)]];
else
{
	$odor = find_odorant(@$_REQUEST['m']);
	if (!$odor)
	{
		header("Location: odorants.php");
		exit;
	}
}

chdir(__DIR__);
chdir("..");
$fullname = $odor['full_name'];
ensure_sdf_exists($fullname);
$sdfname = str_replace(' ','_',"sdf/$fullname.sdf");

$forms = [];
if (@$odor['isomers']) $forms = array_merge($forms, $odor['isomers']);
if (@$odor['forms']) $forms = array_merge($forms, $odor['forms']);
if (count($forms))
{
    if (@$_REQUEST['iso'] || @$_REQUEST['form'])
    {
        $form = @$_REQUEST['iso'] ?: @$_REQUEST['form'];
        if (isset($odor["preiso"]) || isset($odor["preform"]))
        {
            $l = strlen(@$odor["preiso"] ?: @$odor["preform"]);
            $formname = substr($fullname, 0, $l).$form."-".substr($fullname, $l);
            $sdfname = str_replace(' ','_',"sdf/$formname.sdf");
        }
        else $sdfname = str_replace(' ','_',"sdf/$form-$fullname.sdf");
    }
    else if (isset($odor['isomers'])) $sdfname = str_replace(' ','_',"sdf/".(array_keys($odor["isomers"])[0])."-$fullname.sdf");
}
if (file_exists($sdfname))
{
    $sdfdat = file_get_contents($sdfname);
}
else die("File not found $sdfname.\n");

$sdfdat = explode("\n", $sdfdat);
$sdfdat[0] = $fullname;
$sdfdat = implode("\n", $sdfdat);

echo $sdfdat;
