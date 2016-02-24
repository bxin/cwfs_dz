#!/usr/bin/env python

# @author: Bo Xin
# @      Large Synoptic Survey Telescope

import mysql.connector
cnx = mysql.connector.connect(user='root', password='mariadb',
                              host='140.252.32.27',
                              database='WFSCatalog')
cursor = cnx.cursor()

query = ("select * from BrightStars where ra>%s and ra<%s")

cursor.execute(query, ('127.5349', '128.55'))

fidw = open('output/testMysql.txt', 'w')

for (brightStarID, starID, ra, decl, muRA, muDecl, magB,
     magV, magU, magG, magR, magI, magZ, magY, magJ, magH, magK,
     w1, w2, w3, w4, magSST, flag) in cursor:
    fidw.write('%ld, %9.4f, %9.4f\n'%(brightStarID, ra, decl))

cursor.close()
cnx.close()
fidw.close()
