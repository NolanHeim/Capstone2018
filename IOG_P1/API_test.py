from flask import Flask, url_for, request
import uuid
import json

app = Flask(__name__)

@app.route('/')
def api_root():
    return 'Welcome'

@app.route('/articles')
def api_articles():
    
    print("AMBROSIA")
    #window = ["052996 3pm", "052996 4pm"]
    #sat_id = uuid.uuid4()
    target = [49.0, -120.0]

    #the poi json
    poi = {
        "startTime": "052896 3pm",
        "endTime": "052996 10pm",
        }
    
    #the visibility input json
    dictionary = {
        "Target" : target,
        "POI" : poi,            
    } 
    
    response = json.dumps(dictionary)
    return response    
    
    
@app.route('/articles/<articleid>')
def api_article(articleid):
    print("GALA")
    return 'You are reading ' + articleid

@app.route('/visibility/search', methods = ['POST'])
def api_search():
    #print(json.dumps(request.json))
    return json.dumps(request.json)

if __name__ == '__main__':
    app.run()