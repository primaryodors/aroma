<?php

function output_dlmenu_style()
{
?>
#dlmenu
{
    position: absolute;
    box-shadow: 25px 25px 35px rgba(0,0,0,0.5);
    z-index: 10000;
    background: #234;
    padding: 20px;
    font-size: small;
}

#txtview, #htmview
{
    position: absolute;
    z-index: 10000;
    width: 80%;
    min-height: 80%;
    margin: 5% 10%;
    top: 5%;
    border-top: 2px solid #68a;
    border-left: 2px solid #468;
    background-color: #246;
    border-right: 2px solid #123;
    border-bottom: 2px solid #012;
    box-shadow: 0px 0px 50px #000, 0px 0px 100px #000, 0px 0px 200px #000;
}

#htmview
{
    padding: 13px;
    padding-top: 4px;
}

#txttext, #htmhtml
{
    width: 100%;
    height: 100%;
}

.closebtn
{
    display: inline-block;
    border-top: 1px solid #ff0;
    border-left: 1px solid #f90;
    background-color: #f00;
    border-right: 1px solid #900;
    border-bottom: 1px solid #300;
    color: #fff!important;
    text-decoration: none!important;
    cursor: pointer;
    padding: 0px 5px;
}

.ctxmenu
{
    margin: 0px;
    padding-inline-start: 0px;
}

.ctxmenu li
{
    display: block;
    margin-bottom: 2px;
}

.ctxmenu a
{
    display: inline-block;
    border-top: 2px solid #68a;
    border-left: 2px solid #468;
    background-color: #246;
    border-right: 2px solid #123;
    border-bottom: 2px solid #012;
    padding: 0px 5px;
    font-size: small;
}
<?php
}

function output_dlmenu_script()
{
    ?>
function show_dlmenu(e, prot, lig, v, d, abn = "", ibn = "")
{
    var dlmenu = $("#dlmenu")[0];

    $("#dl_mnu_prot")[0].innerText = prot;
    $("#dl_mnu_lig")[0].innerText = decodeURIComponent(lig);
    $("#dl_cd_v")[0].innerText = v;
    if (d == "pd") $("#dl_dock")[0].innerText = "AromaDock";
    else if (d == "vina") $("#dl_dock")[0].innerText = "AutoDock Vina";
    else $("#dl_dock")[0].innerText = "(unknown)";
    $("#dl_acv_mdl")[0].setAttribute("href", "download.php?obj=model&prot="+prot+"&odor="+lig+"&mode=active"+abn);
    $("#vw_acv_mdl_3d")[0].setAttribute("href", "viewer.php?view=pred&prot="+prot+"&odor="+lig+"&mode=active+abn");
    $("#dl_iacv_mdl")[0].setAttribute("href", "download.php?obj=model&prot="+prot+"&odor="+lig+"&mode=inactive"+ibn);
    $("#vw_iacv_mdl_3d")[0].setAttribute("href", "viewer.php?view=pred&prot="+prot+"&odor="+lig+"&mode=inactive"+ibn);
    $("#dl_acv_dc")[0].setAttribute("href", "download.php?obj=dock&prot="+prot+"&odor="+lig+"&mode=active"+abn);
    $("#vw_acv_dc_3d")[0].setAttribute("href", "viewer.php?view=dock&prot="+prot+"&odor="+lig+"&mode=active"+abn);
    $("#dl_iacv_dc")[0].setAttribute("href", "download.php?obj=dock&prot="+prot+"&odor="+lig+"&mode=inactive"+ibn);
    $("#vw_iacv_dc_3d")[0].setAttribute("href", "viewer.php?view=dock&prot="+prot+"&odor="+lig+"&mode=inactive"+ibn);
    $("#dl_json")[0].setAttribute("href", "download.php?obj=json&prot="+prot+"&odor="+lig);

    var x = e.pageX + 5, y = e.pageY;
    if (y > window.innerHeight - 400) y = window.innerHeight - 400;

    dlmenu.style.left = `${x}px`;
    dlmenu.style.top = `${y}px`;

    $(dlmenu).show();
}

function view_file(url)
{
    var txt = $('#txttext')[0];
    var vw = $('#txtview');
    txt.value = "Loading...";

	$.ajax(
	{
		url: url,
		cache: false,
		success: function(result)
		{
            txt.value = result;
            vw.show();
            $('#txttext').css('height', parseInt($('#txttext')[0].parentElement.getClientRects()[0].height));
        }
    });
    vw.show();
}

function view_html_file(url)
{
    var htm = $('#htmhtml')[0];
    var vw = $('#htmview');
    htm.innerHTML = "Loading...";

	$.ajax(
	{
		url: url,
		cache: false,
		success: function(result)
		{
            htm.innerHTML = result;
            vw.show();
        }
    });
    vw.show();
}
<?php
}

function output_dlmenu_div()
{   ?>
<div id="dlmenu" onclick="$('#dlmenu').hide();">
    <div style="margin-bottom: -25px; text-align: right;">
        <a href="#" class="closebtn" onclick="$('#dlmenu').hide();">&#xd7;</a>
    </div>
    <div id="dl_mnu_prot">&nbsp;</div>
    <div id="dl_mnu_lig">&nbsp;</div>
    <br>
    Code version:<div id="dl_cd_v">&nbsp;</div><br>
    Docker used:<div id="dl_dock">&nbsp;</div><br>
    Files:<br>
    <table class="ctxmenu">
        <tr><td>Active model:</td>
            <td><a href="#" onclick="view_file($('#dl_acv_mdl')[0].href); $('#dlmenu').hide();">text</a></td>
            <td><a id="vw_acv_mdl_3d" href="" target="_3d">3D</a></td>
            <td><a id="dl_acv_mdl" href="" target="_dl" onclick="$('#dlmenu').hide();">download</a></td>
        </tr>
        <tr><td>Inactive model:</td>
            <td><a href="#" onclick="view_file($('#dl_iacv_mdl')[0].href); $('#dlmenu').hide();">text</a></td>
            <td><a id="vw_iacv_mdl_3d" href="" target="_3d">3D</a></td>
            <td><a id="dl_iacv_mdl" href="" target="_dl" onclick="$('#dlmenu').hide();">download</a></td>
        </tr>
        <tr><td>Active dock:</td>
            <td><a href="#" onclick="view_html_file($('#dl_acv_dc')[0].href.replace('download.php','tabbed.php')); $('#dlmenu').hide();">text</a></td>
            <td><a id="vw_acv_dc_3d" href="" target="_3d">3D</a></td>
            <td><a id="dl_acv_dc" href="" target="_dl" onclick="$('#dlmenu').hide();">download</a></td>
        </tr>
        <tr><td>Inactive dock:</td>
            <td><a href="#" onclick="view_html_file($('#dl_iacv_dc')[0].href.replace('download.php','tabbed.php')); $('#dlmenu').hide();">text</a></td>
            <td><a id="vw_iacv_dc_3d" href="" target="_3d">3D</a></td>
            <td><a id="dl_iacv_dc" href="" target="_dl" onclick="$('#dlmenu').hide();">download</a></td>
        </tr>
        <tr><td>JSON entry:</td>
            <td><a href="#" onclick="view_file($('#dl_json')[0].href); $('#dlmenu').hide();">text</a></td>
            <td>&nbsp;</td>
            <td><a id="dl_json" href="" target="_dl" onclick="$('#dlmenu').hide();">download</a></td>
        </tr>
    </table>
</div>
<script>
$('#dlmenu').hide();
</script>

<div id="txtview" style="display: none;">
    <div style="text-align: right;">
        <a href="#" class="closebtn" onclick="$('#txtview').hide();">&#xd7;</a>
    </div>
    <textarea id="txttext" enabled="false">
        Loading...
    </textarea>
</div>

<div id="htmview" style="display: none;">
    <div style="text-align: right;">
        <a href="#" class="closebtn" onclick="$('#htmview').hide();">&#xd7;</a>
    </div>
    <div id="htmhtml">
        Loading...
    </div>
</div>

<?php
}