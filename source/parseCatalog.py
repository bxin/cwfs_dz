#!/usr/bin/env python

# @author: Bo Xin
# @      Large Synoptic Survey Telescope

import mysql.connector
cnx = mysql.connector.connect(user='root', password='mariadb',
                              host='140.252.32.27',
                              database='brightStars')
cursor = cnx.cursor()

query = ("select * from BrightStars where ra>%9.4f and ra<%9.4f")

cursor.execute(query, (128.5349, 128.535))

for (brightStarID, starID, ra, decl, muRA, muDecl, magB,
     magV, magU, magG, magR, magI, magZ, magY, magJ, magH, magK,
     w1, w2, w3, w4, magSST, flag) in cursor:
    print('%ld, %9.4f, %9.4f\n'%(brightStarID, ra, decl))

cursor.close()
cnx.close()
