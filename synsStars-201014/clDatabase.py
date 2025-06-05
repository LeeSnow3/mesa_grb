# Written by David Komanek and Richard Wunsch
# History:
#   April 2019 ... initiated as database/db2windcalc
#   July  2019 ... converted to library clDatabase, generalised, cleaned

import MySQLdb
import argparse

#connect to database
def dbConnect(user, passwd):
    try:
        conn = MySQLdb.connect(host='127.0.0.1' \
        , user = user, passwd = passwd, db = 'massive-clusters')
    except:
        print('# Error: unable to connect to the database.')
        print('# Check the username and password.')
        print('# MySQL database should be accessible at localhost:3306.')
        print('# Use e.g.: ssh -N -L 3306:localhost:3306 richard@galaxy.asu.cas.cz')
        conn = None
    return conn

def dbDisconnect(conn):
    conn.cursor().close()
    conn.close()
    return

# retrieves list of cluster names from a specific table
def dbGetGCNames(conn, table):
    cursor = conn.cursor()
    cursor.execute("SELECT Name FROM {}".format(table))
    names = cursor.fetchall()
    return sum(names, ())

# retrieves mass of a cluster "name" from a given table
def dbGetMGC(conn, name, table = 'BaumgardtHilker18'):
    cursor = conn.cursor()
    n = cursor.execute("SELECT Mass FROM {} WHERE Name = \'{}\'".format(table, name))
    mass = cursor.fetchall()
    if len(mass) == 1:
        mass = mass[0][0]
    else:
        mass = None


    # if name not found, print list of available GC and their masses
    if not mass:
        names = dbGetGCNames(conn, 'BaumgardtHilker18')
        print('List of Galactic Globular Clusters with masses in BaumgardtHilker18:')
        for n in names:
            m = dbGetMGC(conn, n)
            print('  {:20} ... {:10.1f} MSun'.format(n, m))

    return mass

def dbGetN1Ntot(conn, name, table = 'Milone16'):
    cursor = conn.cursor()
    n = cursor.execute ("SELECT N1_Ntot FROM {} WHERE Name = \'{}\'".format(table, name))
    N1Ntot = cursor.fetchall()
    if len(N1Ntot) == 1:
        N1Ntot = N1Ntot[0][0]
    else:
        N1Ntot = None

    # if not found, print the list of clusters that have both mass and N1/Ntot
    if not N1Ntot:
        namesBH18 = set(dbGetGCNames(conn, 'BaumgardtHilker18'))
        namesM16  = set(dbGetGCNames(conn, 'Milone16'))
        intersect = namesBH18.intersection(namesM16)
        print('List of GGC in both BaumgardtHilker18 and Milone16:')
        for n in intersect:
            m = dbGetMGC(conn, n)
            r = dbGetN1Ntot(conn, n)
            print('  {:20} ... {:10.1f} MSun ... {: 7.5f}'.format(n, m, r))

    return N1Ntot

def dbGetRGC(conn, name, table = 'BaumgardtHilker18'):
    cursor = conn.cursor()
    n = cursor.execute ("SELECT r_c, r_hm FROM {} WHERE Name = \'{}\'".format(table, name))
    radii = cursor.fetchall()
    if len(radii) == 1:
        Rc  = radii[0][0]
        Rhm = radii[0][1]
    else:
        Rc, Rhm = (None, None)
    return (Rc, Rhm)

# retrieves metallicity of a cluster "name" from a given table
def dbGetZ(conn, name, table = 'Krause16'):
    cursor = conn.cursor()
    n = cursor.execute("SELECT Fe_H FROM {} WHERE Name = \'{}\'".format(table, name))
    Z = cursor.fetchall()
    if len(Z) == 1:
        Z = Z[0][0]
    else:
        Z = None


    # if name not found, print list of available GC and their metallicities
    if not Z:
        names = dbGetGCNames(conn, 'Krause16')
        print('List of Galactic Globular Clusters with metallicities in Krause16:')
        for n in names:
            Z = dbGetZ(conn, n)
            print('  {:30} ... {:5.2f}'.format(n, Z))

    return Z

# retrieves age of a cluster "name" from a given table
def dbGetAge(conn, name, table = 'Krause16'):
    cursor = conn.cursor()
    n = cursor.execute("SELECT Age FROM {} WHERE Name = \'{}\'".format(table, name))
    Age = cursor.fetchall()
    if len(Age) == 1:
        Age = Age[0][0]
    else:
        Age = None


    # if name not found, print list of available GC and their age
    if not Age:
        names = dbGetGCNames(conn, 'Krause16')
        print ('List of Galactic Globular Clusters with age in Krause16:')
        for n in names:
            Age = dbGetAge(conn, n)
            print ('  {:30} ... {:5.2f}'.format(n, Age))

    return Age

# retrieves C5 of a cluster "name" from a given table
def dbGetC5(conn, name, table = 'Krause16'):
    cursor = conn.cursor()
    n = cursor.execute("SELECT C_5 FROM {} WHERE Name = \'{}\'".format(table, name))
    C5 = cursor.fetchall()
    if len(C5) == 1:
        C5 = C5[0][0]
    else:
        C5 = None


    # if name not found, print list of available GC and their C5
    if not C5:
        names = dbGetGCNames(conn, 'Krause16')
        print ('List of Galactic Globular Clusters with age in Krause16:')
        for n in names:
            C5 = dbGetC5(conn, n)
            print ('  {:30} ... {:5.2f}'.format(n, C5))

    return C5

# retrieves retention fraction of a cluster "name" from a given table
def dbGetRet(conn, name, table = 'Baumgardt19'):
    cursor = conn.cursor()
    n = cursor.execute("SELECT Mass, InitialMass FROM {} WHERE Name = \'{}\'".format(table, name))
    masses = cursor.fetchall()
    if len(masses) == 1:
        ret = masses[0][0] / masses[0][1]
    else:
        ret = None


    # if name not found, print list of available GC and their C5
    if not ret:
        names = dbGetGCNames(conn, table)
        print ('List of Galactic Globular Clusters with retention fraction in Krause16:')
        for n in names:
            ret = dbGetRet(conn, n)
            print ('  {:30} ... {:5.2f}'.format(n, ret))

    return ret

def getParser():
# {{{
    """Parses arguments using argparse module."""

    # setup parser
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--dbuser", default="", type=str, help="Username for retrieving an object from the database. \
    Default is ''. Database 'massive-clusters' must be accessible on localhost:3306")
    parser.add_argument("--dbpass", default="", type=str, help="Password for retrieving an object from the database. Default is ''.")
    parser.add_argument("--dbname", default="", type=str, help="Name of the object whose mass \
    will be obtained from the database. Default is ''. If the object is not \
    found, it prints list of all clusters in the BaumgardtHilker18 table")
    return parser
# }}}
