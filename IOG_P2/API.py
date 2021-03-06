# -*- coding: utf-8 -*-
#
# Api.py
#
# (PHASE 2)
#
# This file contains all the code needed for the REST API that wraps the project.
# It instantiates the Mission Creator (top level file for the program) and can be
# communicated with via the command line. See the final report Appendix for full
# running instructions.
#
# Initial Creation Date: 11/03/2018
#
# Written by Jordan Jones and Nolan Heim
#

from flask import Flask, request, Response
import json
import uuid
import os
from Mission_Creator_REST import *

api = Flask(__name__)


#Base
@api.route('/')
def api_root():
    return "Welcome to the Imaging Opportunity Generator"


#Post
@api.route('/opportunity/search', methods = ['POST'])
def api_search():
        
    request_uuid = uuid.uuid4()
    status = 500

    opportunities = MC.generate_imaging_opportunities(request.json, request_uuid)        
      
    if opportunities == "ERR": #Mission Creator failed due to invalid params
        response = "ERROR - Wrong/missing input object attribute(s)"
        #Create the Response object
        resp = Response(response, status=400, mimetype='application/json')
        return resp
    
    output_json = {
        "id":str(request_uuid),
        "Opportunities":opportunities
    }

    # Saves output for use by GET function        
    with open(results_path + str(request_uuid) + ".json", "w") as save_output:
        save_output.write(json.dumps(output_json))
        status = 200
    
    response = json.dumps(output_json)

    #Create the Response object    
    resp = Response(response, status=status, mimetype='application/json')

    return resp


#Get
@api.route('/opportunity/<articleid>', methods = ['GET'])
def api_get_results(articleid):

    #no such results id exists, return code 400
    if not (os.path.isfile(results_path+str(articleid)+".json")):
        response = "ERROR - No results for specified ID"
        #Create the Response object
        resp = Response(response, status=400)
        return resp
    
    else: #results for given ID do exist
        with open(results_path+str(articleid)+".json", "r") as read_file:
            response = read_file.readlines()
        
        #Create the Response object
        resp = Response(response, status=200, mimetype='application/json')
        return resp
    
    response = "ERROR - Backend error"
    #Create the Response object
    resp = Response(response, status=500)
    return resp


# Code to actually run the API. Then open a separate command window to communicate with it.
if __name__ == '__main__':
    datapath = "../../Data/"
    results_path = "../../Saved Results/"
    parsed_datapath = "../../Parsed Data/"
    MC = Mission_Creator_REST(parsed_datapath, results_path)
    api.run()