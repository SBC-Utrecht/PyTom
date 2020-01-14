def messageResponse(message,code = 400):
    """
    messageResponse: Generates a html message on the fly
    """

def htmlHeader200(server):
    server.send_response(200)
    server.send_header('Host','PyTom Frontend Server')
    server.send_header("Content-type","text/html")
    server.end_headers()
    
    
def pngHeader200(server):
    server.send_response(200)
    server.send_header('Host','PyTom Frontend Server')
    server.send_header("Content-type","image/png")
    server.end_headers()