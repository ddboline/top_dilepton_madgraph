#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import glob

def process_files( arg , dirname , names ) :
    for name in names :
        if 'COMPONENTS_F' in name :
            #print dirname , name
            for line in open( '%s/%s' % ( dirname , name ) ) :
                os.system( 'gfortran -fpic -g -ffixed-line-length-132 -c %s/%s.f -o obj/%s.o' % ( dirname , line.strip() , line.strip() ) )

os.system( 'mkdir -p obj lib' )
os.system( 'rm obj/*.o lib/*.so' )
os.path.walk( '.' , process_files , None )
os.system( 'gfortran -shared -o lib/libtop_dilepton_madgraph.so obj/*.o' )