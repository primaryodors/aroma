<?php

function get_3d_distance($pt1, $pt2)
{
    return sqrt( pow($pt1[0]-$pt2[0], 2) + pow($pt1[1]-$pt2[1], 2) + pow($pt1[2]-$pt2[2], 2) );
}

function get_2d_distance($pt1, $pt2)
{
    return sqrt( pow($pt1[0]-$pt2[0], 2) + pow($pt1[1]-$pt2[1], 2) );
}

function scale_point($pt, $scale)
{
    $r = sqrt( pow($pt[0], 2) + pow($pt[1], 2) + pow($pt[2], 2) );
    if ($r)
    {
        return [ $pt[0]*$scale/$r, $pt[1]*$scale/$r, $pt[2]*$scale/$r ];
    }
    else return [0,0,0];
}
function rotate3D($point, $source, $axis, $theta)
{
    // Originally from https://web.archive.org/web/20131229124319/http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/

    $axisr = sqrt( pow($axis[0], 2) + pow($axis[1], 2) + pow($axis[2], 2) );
    if (!$axisr) return $point;
    $lvec = $axis;

    list($x,$y,$z) = $point;
    $u = $lvec[0] / $axisr;
    $v = $lvec[1] / $axisr;
    $w = $lvec[2] / $axisr;

    list($a,$b,$c) = $source;

    $u2 = $u*$u;
    $v2 = $v*$v;
    $w2 = $w*$w;
    $sint = sin($theta);
    $cost = cos($theta); 
    $_1_cost = (1.0 - $cost);

    $x1 = ($a * ($v2+$w2) - $u * ($b*$v + $c*$w - $u*$x - $v*$y - $w*$z)) * $_1_cost
               + $x * $cost
               + (-$c*$v + $b*$w - $w*$y + $v*$z) * $sint;

    $y1 = ($b * ($u2+$w2) - $v * ($a*$u + $c*$w - $u*$x - $v*$y - $w*$z)) * $_1_cost
               + $y * $cost
               + ( $c*$u - $a*$w + $w*$x - $u*$z) * $sint;

    $z1 = ($c * ($u2+$v2) - $w * ($a*$u + $b*$v - $u*$x - $v*$y - $w*$z)) * $_1_cost
               + $z * $cost
               + (-$b*$u + $a*$v - $v*$x + $u*$y) * $sint;

    return [$x1,$y1,$z1];
}

function find_3d_angle($A, $B, $source)
{
    $lA = scale_point([ $A[0]-$source[0], $A[1]-$source[1], $A[1]-$source[2] ], 1);
    $lB = scale_point([ $B[0]-$source[0], $B[1]-$source[1], $B[2]-$source[2] ], 1);

    // https://stackoverflow.com/questions/1211212/how-to-calculate-an-angle-from-three-points
    $P12 = get_3d_distance($lA, [0,0,0]);
    $P13 = get_3d_distance($lB, [0,0,0]);
    $P23 = get_3d_distance($lA, $lB);

    $param = ($P12*$P12 + $P13*$P13 - $P23*$P23)/(2 * $P12 * $P13+.00000000001);
    if ($param < -1) $param = -1;
    if ($param >  1) $param =  1;
    $retval = acos($param);
    if (is_nan($retval))
    {
        echo "P12 $P12 P13 $P13 P23 $P23\n";
        exit;
    }
    return $retval;
}

function compute_normal($pt1, $pt2, $pt3)
{
    $U = [ $pt2[0] - $pt1[0], $pt2[1] - $pt1[1], $pt2[2] - $pt1[2] ];
    $V = [ $pt3[0] - $pt1[0], $pt3[1] - $pt1[1], $pt3[2] - $pt1[2] ];

    return
    [	$U[1] * $V[2] - $U[2] * $V[1],
        $U[2] * $V[0] - $U[0] * $V[2],
        $U[0] * $V[1] - $U[1] * $V[0]
    ];
}

function align_points_3d($point, $align, $center)
{
    $n = compute_normal($point, $align, $center);
    $nr = sqrt( pow($n[0], 2) + pow($n[1], 2) + pow($n[2], 2) );

    if ($nr < 0.0001)
    {
        $lpt = scale_point($point, 1);
        $lan = scale_point($align, 1);

        $result = [];
        if (get_3d_distance($lpt, $lan) < 0.01)
        {
            $result = $n;
            $result[3] = 0;
            return $result;
        }

        $pt = [0,0,1];
        $n = compute_normal($point, $align, $pt);
        $nr = sqrt( pow($n[0], 2) + pow($n[1], 2) + pow($n[2], 2) );
        if ($nr < 0.1)
        {
            $pt = [0,0,1];
            $n = compute_normal($point, $align, $pt);
        }

        $result = $n;
        $result[3] = M_PI;
        return $result;
    }

    // Find the 3D angle between pp and pl relative to center.
    $theta = find_3d_angle($point, $align, $center);
    // cout << " theta = " << theta << " ";

    // Rotate pl positively or negatively along that normal by the found angle, and use the better of the two values.
    $plus  = rotate3D($point, $center, $n,  $theta);
    $minus = rotate3D($point, $center, $n, -$theta);

    $rplus  = get_3d_distance($plus, $align);
    $rminus = get_3d_distance($minus, $align);

    if ($rplus <= $rminus) $angle =  $theta;
    else                   $angle = -$theta;

    $rot = $n;
    $rot[3] = $angle;
    return $rot;
}

