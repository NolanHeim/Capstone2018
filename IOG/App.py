# -*- coding: utf-8 -*-
#
# App.py
#
# Initial Creation Date: 11/03/2018
#
# Written by Jordan Jones and Nolan Heim
#

from flask import Flask, url_for, request
import json
import uuid
import os
from MissionCreatorRest import *

app = Flask(__name__)


def __init__(self):
    self.datapath = "../../Data/"
    self.basepath = "../../SavedResults/"
    self.parsed_datapath = "../../Parsed Data/"
    self.MC = MissionCreatorREST(self.datapath)

@app.route('/')
def api_root():
    return 'Welcome'




#post
@app.route('/visibility/search', methods = ['POST'])
def api_search():
        
    request_uuid = uuid.uuid4()
    opportunities = self.MC.generate_imaging_opportunities(request.json, request_uuid)        

        
    output_json = {
        'id':request_uuid,
        'Opportunities':opportunities
    }
        
    response = json.dump(output_json)
    response.status_code = 200
        
    return response


#get
@app.route('/visibility/<id>', methods = ['GET'])
def api_getresults(articleid):
    print('You are reading ' + articleid)
    
    #no such results id exists
    if not (os.path.isfile(self.basepath + articleid + ".json")):
        response = ""
        response.status_code = 400
        return response
    
    with open(self.basepath + articleid + ".json", "r") as read_file:
        response = read_file
        response.status_code = 200
        return response
    
    response = ""
    response.status_code = 500
    return response

if __name__ == '__main__':
    app.run()