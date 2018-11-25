# -*- coding: utf-8 -*-
#
# App.py
#
# Initial Creation Date: 11/03/2018
#
# Written by Jordan Jones and Nolan Heim
#

from flask import Flask, url_for, request, Response
import json
import uuid
import os
from Mission_Creator_REST import *

api = Flask(__name__)


@api.route('/')
def api_root():
    return "Welcome to the Imaging Opportunity Generator"


@api.route('/hello')
def api_hello():
    return "Hello"


#post
@api.route('/visibility/search', methods = ['POST'])
def api_search():
        
    request_uuid = uuid.uuid4()
    status = 500

    opportunities = MC.generate_imaging_opportunities(request.json, request_uuid)        
      
    if opportunities == "ERR":
        response = "ERROR - Wrong/missing input object attribute(s)"
        resp = Response(response, status=400, mimetype='application/json')
        return resp
    
    output_json = {
        "id":str(request_uuid),
        "Opportunities":opportunities
    }

    #saves output for use by GET function        
    with open(results_path + str(request_uuid) + ".json", "w") as save_output:
        save_output.write(json.dumps(output_json))
        status = 200
    
    response = json.dumps(output_json)
    
    resp = Response(response, status=status, mimetype='application/json')

    return resp


#get
@api.route('/visibility/<articleid>', methods = ['GET'])
def api_get_results(articleid):

    #no such results id exists
    if not (os.path.isfile(results_path+str(articleid)+".json")):
        response = "ERROR - No results for specified ID"
        resp = Response(response, status=400)
        return resp
    
    else:
        with open(results_path+str(articleid)+".json", "r") as read_file:
            response = read_file.readlines()

        resp = Response(response, status=200, mimetype='application/json')
        return resp
    
    response = "ERROR - Backend error"
    resp = Response(response, status=500)
    return resp


if __name__ == '__main__':
    datapath = "../../Data/"
    results_path = "../../Saved Results/"
    parsed_datapath = "../../Parsed Data/"
    MC = Mission_Creator_REST(parsed_datapath, results_path)
    parser = Parser(datapath, "")
    #parser.parse_data(parsed_datapath)    
    api.run()