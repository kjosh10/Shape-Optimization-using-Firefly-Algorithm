

# In[83]:
"""
importing necessary dependencies
"""
import numpy as np
import matplotlib.pyplot as plt

def pointCalc(point_array, t, cubic=False):
    """
    method to compute point using Bezier Curve

    :param point_array: list, a list of points
    :param t: float, a parameter for Bezier Curve
    :param cubic: boolean, a boolean for using either Quadratic or
                    cubic Bezier curve
    :return points_x: flaot, x co-ordinate of point on Bezier curve
    :return points_y: flaot, y co-ordinate of point on Bezier curve
    """
    if (cubic):
        points_x = point_array[0][0]*(-t**3+3*t**2-3*t+1) +point_array[1][0]*(3*t**3-6*t**2+3*t)+point_array[2][0]*(-3*t**3+3*t**2)+point_array[3][0]*t**3
        points_y = point_array[0][1]*(-t**3+3*t**2-3*t+1)+point_array[1][1]*(3*t**3-6*t**2+3*t)+point_array[2][1]*(-3*t**3+3*t**2)+point_array[3][1]*t**3
        print(points_x, points_y)
    else:                
        points_x = point_array[0][0]*(1-t)**2 +point_array[1][0]*(1-t)*2*t + point_array[2][0]*(t)**2
        points_y = point_array[0][1]*(1-t)**2 +point_array[1][1]*(1-t)*2*t + point_array[2][1]*(t)**2
    return points_x, points_y

def bezierCurveTopPoint(y__middle_curve, point_array):
    """
    method to compute the top point of the Bezier Curve

    :param y__middle_curve: float, y co-ordinate of the middle point
    :return float, top point of the Bezier Curve
    """
    t = 0.5
    return (y__middle_curve - point_array[0][1]*(1-t)**2 - 
            point_array[1][1]**(t)**2)/((1-t)*2*t)

