# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 18:37:34 2016
@author: fernandi
"""


import sys
import shutil
import math
import os

if __name__ == "__main__":
    # args: base_props_file number_robots init_x init_y delta_position
    if len(sys.argv) != 6:
        sys.exit("Error: wrong number of arguments\n" +
                 "Usage: python positionRobots base_props_file "
                 "number_robots init_x init_y delta_position")
    path_config = "../config/"
    dest = sys.argv[1].split(".properties")[0] + "_with_positions.properties"
    if os.path.exists(path_config + dest):
        os.remove(path_config + dest)
    shutil.copyfile(path_config + sys.argv[1], path_config + dest)

    n = int(sys.argv[2])  # number of robots
    init_x = int(sys.argv[3])  # initial x position
    init_y = int(sys.argv[4])  # initial y position
    delta = int(sys.argv[5])  # delta (both x and y) between robots

    robots_per_side = math.ceil(math.sqrt(n))  # square disposition

    print dest, n, init_x, init_y, delta

    x = init_x
    y = init_y
    count_r = 0
    orient = 90
    s = "\n\n\n\n################\n#Robot positions\n#############\n"

    with open(path_config + dest, "a") as copied_file:
        for i in xrange(n):
            s += "robot[" + str(i) + "].x = " + str(x) +\
                "\nrobot[" + str(i) + "].y = " + str(y) + "\n" +\
                "robot[" + str(i) + "].orientation = " + str(orient) +\
                "\t# 0...359, clockwise -- default is 0\n" +\
                "robot[" + str(i) + "].groupId=0 \t# default is 0\n"
            count_r += 1
            if count_r == robots_per_side:
                count_r = 0
                x = init_x
                y += delta
            else:
                x += delta
        print s
        copied_file.write(s)

    # read arg pos 1er
    # read arg delta (distance inter-robots)
    # copy file properties

    # get nb robots
    # for each robot set initial position
    # set random orientation?
