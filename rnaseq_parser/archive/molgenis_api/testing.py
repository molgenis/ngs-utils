# -*- coding: utf-8 -*-
import json
import cProfile
from molgenis_api import security
from molgenis_api import molgenis
import uuid
import pstats
import random
from io import StringIO

import requests
from datetime import datetime


def test_duplicateValues():
    with molgenis.Connect_Molgenis('http://localhost:8080') as connection:
        print(connection.add_entity_rows('test_duplicateTestTable',[{'id':5,'info':'duplicate row?'},{'id':6,'info':'not duplicate'}]))











test_duplicateValues()








def profile_multirows():
    with molgenis.Connect_Molgenis('http://localhost:8080') as connection:
        connection.delete_all_entity_rows('test_testTable1')
        connection.delete_all_entity_rows('test_testTable2')
        connection.delete_all_entity_rows('test_testTable3')
        #to_add = []
        pr = cProfile.Profile()
        pr.enable()
        #for i in range(0,1000):
        #    connection.add_entity_rows('test_testTable2',[{'example':'row_'+str(i)}])
        #pr.disable()
        #s = StringIO()
        #sortby = 'cumulative'
        #ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        #ps.print_stats()
        #print (s.getvalue())
        #pr.enable()
        to_add = []
        for i in range(0,1000):
            to_add.append({'example':'row_'+str(i)})
        connection.add_entity_rows('test_testTable2',to_add)
        pr.disable()
        s = StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print (s.getvalue())

def old_profiling():
    connection = molgenis.Connect_Molgenis('https://molgenis39.target.rug.nl').__enter__()
    #connection.delete_all_entity_rows('public_rnaseq_2_Samples')
    connection.delete_all_entity_rows('public_rnaseq_2_Analysis_info')
    connection._added_by_default = True
    connection._add_datetime_default = True
    request_url = 'https://molgenis39.target.rug.nl/api/v1/public_rnaseq_2_Analysis_info'
    
    def post_10_api():
        for i in range(0,10):
            connection.add_entity_row('public_rnaseq_2_Analysis_info',data={'analysis_description':'profiling api','id':random.randint(2, 9999999)})
    def get_10_api():
        for i in range(0,10):
            connection.get_entity('public_rnaseq_2_Analysis_info')
    
    def post_10():
        for i in range(0,10):
            requests.post(request_url, data=json.dumps({'analysis_description':'profiling direct request','id':random.randint(2,9999999),'added_by':True,'datetime_added':str(datetime.now())}), headers=connection.headers)
    
    def get_10():
        for i in range(0,10):
            requests.get(request_url, headers=connection.headers)
    
    functions = ['post_10', 'get_10','post_10_api','get_10_api']
    for f in functions:
        print (f)
        cProfile.run(f+'()', 'restats')
        p = pstats.Stats('restats')
        p.strip_dirs().sort_stats('cumulative').print_stats(25)
        print ('-'*20)
