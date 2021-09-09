#!/usr/bin/env pytom
'''
Created on Jul 7, 2011

@author: Thomas Hrabe
'''


import http.server
from pytom.tools.files import getPytomPath

pytomPath = getPytomPath()

verbose = False
#verbose = True

#define where files are stored according to type
fileTypePaths = {'':'','index.html':'','html':'pages','xsl':'xslt','js':'javascript','png':'images','css':'css','py':'../serverpages'}

def dispatchRequest(server,filename,parameters,verbose=False):
    """
    dispatchRequest: Dispatches a filename to a function
    @param filename:
    @param parameters:  
    """
    
    if verbose:
        print("Dispatching serverpage request : ", filename, parameters)
    
    if filename == 'loadParticleList.py':
        from pytom.frontend.serverpages.loadParticleList import run
    elif filename == 'apiDocumentationURL.py':
        from pytom.frontend.serverpages.apiDocumentationURL import run
    elif filename == 'imageSlice.py':
        from pytom.frontend.serverpages.imageSlice import run
    elif filename == 'createAlignmentJob.py':
        from pytom.frontend.serverpages.createAlignmentJob import run
    elif filename == 'loadAlignmentJob.py':
        from pytom.frontend.serverpages.loadAlignmentJob import run
    elif filename == 'createLocalizationJob.py':
        from pytom.frontend.serverpages.createLocalizationJob import run
    elif filename == 'createMCOEXMXJob.py':
        from pytom.frontend.serverpages.createMCOEXMXJob import run
    elif filename == 'createMCOACJob.py':
        from pytom.frontend.serverpages.createMCOACJob import run
    elif filename == 'loadProjectionList.py':
        from pytom.frontend.serverpages.loadProjectionList import run
    elif filename == 'createReconstructionJob.py':
        from pytom.frontend.serverpages.createReconstructionJob import run
    elif filename == 'sendByteArray.py':
        from pytom.frontend.serverpages.sendByteArray import run
    else:
        return None
    
    
    resultString = run(parameters,verbose)
    return server.sendString(resultString)

def parseQuery(query):
    """
    parseQuery: Will return pairs of parameters seperated by & appended to url starting after the ? . (Common url parsing practice) 
    @author: Thomas Hrabe
    """
    
    from urllib.parse import urlsplit
    
    query = query.replace('%2F','/')
    query = query.replace('%7E','~')
    
    parts = urlsplit(query)
    
    requestParameters = parts[3]
    
    return requestParameters.split('&')

def parsePythonRequest(server):
    """
    parsePythonRequest: Parses a typical http <filename>.py?a=1&b=2 request and sends a string back
    """

    request = server.path
    
    if request.find('?') != -1: 
        splitRequest = request.split('?')
    
        url = splitRequest[0]
        parameters = splitRequest[1]

    else:
        url = request
        parameters = ''
     
    splitURL = url.split('/')
    
    if verbose:
        print(parameters)
        
    fileObject = splitURL[len(splitURL)-1]

    dispatchRequest(server,fileObject,parameters,verbose)

    
    
class PyTomHTTPServer(http.server.BaseHTTPRequestHandler):
    """
    PyTomHTTPServer:
    """
    
    def _returnFilePath(self,path,extension = ''):
        from pytom.tools.files import getPytomPath, checkFileExists
        pytomPath = getPytomPath()
        filePath = pytomPath + '/frontend/html/'
        print(filePath + fileTypePaths[extension] + path)
        if checkFileExists(path):
            return path
        elif checkFileExists(filePath + path):
            return filePath + path
        elif checkFileExists(filePath + fileTypePaths[extension] + path):
            return filePath + fileTypePaths[extension] + path
        elif path in ['','/','/index.html']:
            return filePath + 'index.html'
        else:
            raise IOError('File not found: ' + str(path))
        
        
    def do_ERROR(self):
        """
        do_ERROR: Send error message in case of some wrongdoing
        """
        
    def do_HEAD(self,path,extension):
        """
        do_HEAD
        @param path: URL between / and ? - like index.html or doSomething.py
        @param extension: File extension
        """
        import mimetypes
        from pytom.tools.files import getPytomPath, checkFileExists
        
        if path in ['' ,'/' ,'/index.html']:
            self.send_response(200)
            self.send_header('Host','PyTom Frontend Server')
            self.send_header("Content-type","text/html")
            self.end_headers()
            return
    
        if not extension in fileTypePaths:
            #requested filetype not supported. abort!
            self.send_response(406)
            self.end_headers()
            return
        
        try:
            self._returnFilePath(path, extension)
        except IOError:
            #file not found! abort!
            self.send_response(404)
            self.end_headers()
            return
        
        #respond all is well!
        self.send_response(200)
        self.send_header('Host','PyTom Frontend Server')
        self.send_header("Content-type",mimetypes.guess_type(path)[0])
        self.end_headers()

    def do_GET(self):
        """
        do_GET
        """
        
        from pytom.tools.files import getPytomPath , checkFileExists
        import urllib.parse
        if verbose:
            print('HTTP Request : ',self.path)

        
        parsedURL = urllib.parse.urlparse(self.path)
        pathPartition = parsedURL.path.partition('.') 
        extension = pathPartition[2]
        
        self.do_HEAD(parsedURL.path,extension)
        
        filePath = self._returnFilePath(parsedURL.path, extension)
        
        if extension in ["jpg","png","gif" ,"tiff","ico",".jpg",".png",".gif" ,".tiff",".ico"]:
            self.sendBinary(filePath)
        elif extension == "py":
            parsePythonRequest(self)
        else:
            self.sendStringFile(filePath)
    
    def sendString(self,string):
        self.wfile.write(string)
    
    def sendStringFile(self,filepath):
        
        from pytom.tools.files import readStringFile
        htmlString = readStringFile(filepath)
        self.sendString(htmlString)
        
    def sendBinary(self,filepath):
        
        f = open(filepath,"r")    
        self.wfile.write(f.read())
        f.close()
     
if __name__ == '__main__':
    
    import sys, getopt
    
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.tools.files import checkFileExists
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Start PyTom User Interface.',
                          authors='Thomas Hrabe',
                          options=[ScriptOption(['-h','--hostname'], 'Hostname of server.(name of computer you are logged into)', True, True),
                                   ScriptOption(['-p','--port'], 'Specify port of server (if you do not know what that means, use 8080)', True, True),
                                   ScriptOption(['-b','--browser'], 'Name of browser to start.', True, True),
                                   ScriptOption(['-d','--disableBrowser'], 'Disable browser at start.', False, True),
                                   ScriptOption(['-h', '--help'], 'Help.', False, True)])
 
    defaultHostname = None
    defaultPort = 8080
    defaultBrowser = 'firefox'
    
    if checkFileExists('/Applications/Firefox.app/Contents/MacOS/firefox'):
        defaultBrowser = '/Applications/Firefox.app/Contents/MacOS/firefox'
    
    defaultDisableBrowser = False
    
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        hostname, port, browser, disableBrowser, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()
    if help is True:
        print(helper)
        sys.exit()
    
    if not hostname:
        hostname = defaultHostname
        
    if not port:
        port = defaultPort
    else:
        port = int(port)
        
    if not browser:
        browser = defaultBrowser
    
    if not disableBrowser:
        disableBrowser = defaultDisableBrowser
    
    if not port or not hostname:
        print(helper)
        sys.exit()
        
    print('----------------------------------------------------------------------------------------')        
    print('Starting Pytom Webserver on hostname ' + hostname + ' and port ' + str(port))
    print('Paste the following line into your browser to start')
    print('http://' + hostname + ':' + str(port))
    if not disableBrowser:
        print('Starting browser: ', browser)
    print('')
    print('----------------------------------------------------------------------------------------')
    
    
    httpd = http.server.HTTPServer((hostname,port),PyTomHTTPServer)
    
    
    if not disableBrowser:
        import os
        pid = os.fork()
        
        if pid == 0:
            try:
                print(browser + ' http://' + hostname + ":" + str(port) + '/index.html')
                os.system(browser + ' http://' + hostname + ":" + str(port) )
            except:
                print('Problems starting your browser. Please start browser manually.')
                pass
            sys.exit()
     
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    
    httpd.server_close()
