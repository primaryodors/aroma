
import math
import data.globals

def subtract_points(pt1, pt2):
    result = [0,0,0]
    for i in range(3):
        result[i] = pt1[i] - pt2[i]
    return result

def pt_magnitude(pt):
    return pow(pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2], 0.5)

def scale_point(pt, new_magn):
    m = pt_magnitude(pt)
    if not m: return [0,0,0]
    return [pt[0] / m, pt[1] / m, pt[2] / m]

def pt_distance(pt1, pt2):
    return pt_magnitude(subtract_points(pt1, pt2))

def compute_normal(pt1, pt2, pt3):
    U = subtract_points(pt2, pt1)
    V = subtract_points(pt3, pt1)

    return [ U[1] * V[2] - U[2] * V[1],
             U[2] * V[0] - U[0] * V[2],
             U[0] * V[1] - U[1] * V[0]
           ]

def find_3d_angle(pt1, pt2, center):
    lA = scale_point(subtract_points(pt1, center), 1)
    lB = scale_point(subtract_points(pt2, center), 1)

    P12 = pt_magnitude(lA)
    P13 = pt_magnitude(lB)
    P23 = pt_distance(lA, lB)

    param = (P12*P12 + P13*P13 - P23*P23)/(2 * P12 * P13+.00000000001)
    if param < -1: param = -1
    if param >  1: param =  1

    return math.acos(param)

def rotate3D(pt, center, axis, theta):
    if not pt_magnitude(axis): return pt
    x = pt[0]
    y = pt[1]
    z = pt[2]
    r = pt_magnitude(axis)
    u = axis[0] / r
    v = axis[1] / r
    w = axis[2] / r
    a = center[0]
    b = center[1]
    c = center[2]

    u2 = u*u
    v2 = v*v
    w2 = w*w
    sint = math.sin(theta)
    cost = math.cos(theta)
    cost1 = (1.0 - cost)

    x1 = (a * (v2+w2) - u * (b*v + c*w - u*x - v*y - w*z)) * cost1 \
               + x * cost \
               + (-c*v + b*w - w*y + v*z) * sint

    y1 = (b * (u2+w2) - v * (a*u + c*w - u*x - v*y - w*z)) * cost1 \
               + y * cost \
               + ( c*u - a*w + w*x - u*z) * sint

    z1 = (c * (u2+v2) - w * (a*u + b*v - u*x - v*y - w*z)) * cost1 \
               + z * cost \
               + (-b*u + a*v - v*x + u*y) * sint

    return [x1,y1,z1]

def align_points_3d(point, align, center):
    n = compute_normal(point, align, center)

    if pt_magnitude(n) < 0.000001:
        lpt = scale_point(point)
        lan = scale_point(align)

        if pt_distance(lpt, lan) < 0.0000001:
            return [0, 0, 0, 0]
        else:
            return [n[0], n[1], n[2], math.pi]

    theta = find_3d_angle(point, align, center)

    plus  = rotate3D(point, center, n,  theta)
    minus = rotate3D(point, center, n, -theta)

    rplus  = plus.get_3d_distance(align)
    rminus = minus.get_3d_distance(align)

    if rplus <= rminus: angle =  theta
    else:               angle = -theta

    return [n[0], n[1], n[2], angle]