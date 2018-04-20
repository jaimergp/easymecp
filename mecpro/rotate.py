#!/usr/bin/python

#Copyright (c) 2015 Brigham Young University

#See the file license.txt for copying permission.
        
import math
import operator

def vector(a, b):
    return map(lambda (x,y): y-x, zip(a, b))
    
def add_vec(a, b):
    return map(lambda (x,y): x+y, zip(a, b))

def cross_product(a, b):
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]
    return c
    
def dot_product(a,b):
    return [a[0] * b[0],
            a[1] * b[1],
            a[2] * b[2]]
            
def length(vec):
    return reduce(math.hypot, vec)
    
def unit_vec(vec):
    length = reduce(math.hypot, vec)
    return map(lambda x: x/length, vec)
    
'''Get the RHS rotation axis from three points.
specify the points in the same order/direction you want them rotated'''
def get_rot_axis(a, b, c):
    axis = cross_product(vector(b,a), vector(b,c))
    #normalize the axis
    return unit_vec(axis)
    
def get_rot_angle(a, b):
    a = unit_vec(a)
    b = unit_vec(b)
    return math.acos(length(dot_product(a, b)))
    
def get_rot_matrix(angle, axis):
    c = math.cos(angle)
    t = 1 - c
    s = math.sin(angle)
    ax, ay, az = unit_vec(axis)
    return [      # Rotation Matrix
        [ t*ax*ax + c,     t*ax*ay + s*az,   t*ax*az - s*ay,   0 ],
        [ t*ay*ax - s*az,  t*ay*ay + c,      t*ay*az + s*ax,   0 ],
        [ t*az*ax + s*ay,  t*az*ay - s*ax,   t*az*az + c,      0 ],
        [      0,                 0,                 0,        1 ] ]

def mul_mat_vec(mat, vec):
    result = []
    for row in mat:
        result.append( map(lambda (x,y): x*y, zip(row, vec)) ) # multiply across
    result = map(lambda l: reduce(operator.add, l), result) # add stuff together
    return result
    
'''Rotate a point about an arbitrary axis and point (RHS)
Not sure if this is accurate... But the math at least looks legit.'''
def rotate(point, angle, center, axis):
    rot_matrix = get_rot_matrix(angle, axis)
    vec = vector(center, point)     # rotate the point about the center
    vec.append(0)           # make the vector 4-dimensional
    result = mul_mat_vec(rot_matrix, vec)[:3]     # remove the 4th dimension
    #print point, '-->', add_vec(center, result)
    return add_vec(center, result)  # translate back to the original origin
    
def rotate_many(points, angle, center, axis):
    rot_matrix = get_rot_matrix(angle, axis)
    result = []
    for point in points:
        vec = vector(center, point)     # rotate the point about the center
        vec.append(0)           # make the vector 4-dimensional
        new_point = mul_mat_vec(rot_matrix, vec)[:3]     # remove the 4th dimension
        #print point, '-->', add_vec(center, new_point)
        result.append(add_vec(center, new_point))  # translate back to the original origin
    return result
        
    
if __name__ == "__main__":
    axis = get_rot_axis((1,0,1), (0,1,0), (-1,0,-1))
    print rotate((4,2,1), math.pi, (0,0,0), axis)
