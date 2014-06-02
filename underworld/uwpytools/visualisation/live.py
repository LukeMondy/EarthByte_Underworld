import uwpytools as _uwpytools



def StartVisualisationWebServer(host='127.0.0.1', port=-99999, certfile=None):
    """
    Starts a visualisation web server (via root processor) in a new thread.
    An Init function first needs to be called so that MPI_Init has delegated processor Ids. 

    Args:
        host     (string)        : (optional) The server host address.  Defaults to localhost (127.0.0.1).
        port     (   int)        : (optional) The server listening port.  Defaults to 8999 and autoincrementing if necessary.
        certfile (string)       : (optional) A certificate file (including path) for running an HTTPS web server.

    Returns:
        Nothing
    """
    
    _data = _uwpytools.data()

    if _data == None:
        print "StGermain has not been initialised yet. You need to run one of the Init functions."
        return

    if _uwpytools.rank() == 0:
        defaultport = 8999
        # lets first halt any running instances
        StopVisualisationWebServer()

        if port < 0:
            _port = defaultport
        else:
            _port = port

        import SimpleHTTPServer
        import SocketServer
        import threading
        global _httpd

        SocketServer.TCPServer.allow_reuse_address = True
        Handler = SimpleHTTPServer.SimpleHTTPRequestHandler

        if port < 0:
            trying = True
            while(trying):
                try:
                    _httpd = SocketServer.TCPServer((host, _port), Handler)
                    trying = False
                except Exception, e:
                    _port +=1
                    if _port > defaultport + 10:
                        print "Error. Tried multiple ports without success. Perhaps try setting a port explicitly"
                        print "using the port= option, or shutdown any already running servers."
                        print str(e)
                        raise
        else:
            try:
                _httpd = SocketServer.TCPServer((host, _port), Handler)
            except Exception, e:
                #print "Error. Unable to open port ", _port
                print "Port already being used, continuing as if nothing is wrong. ", _port
                #print str(e)
                #raise
                return

        if certfile:
            import ssl
            _httpd.socket = ssl.wrap_socket(_httpd.socket, certfile=certfile, server_side=True)
        httpd_thread = threading.Thread(target=_httpd.serve_forever)
        httpd_thread.setDaemon(True)
        httpd_thread.start()
        print("Serving visualisation web server at %s:%i" % (host,_port))

        return host+str(port)


def StopVisualisationWebServer():
    """
    Halts any running visualisation servers

    Args:
        None

    Returns:
        Nothing
    """
    global _httpd
    if _httpd:
        _httpd.shutdown()
        _httpd = None
        print "Visualisation web server halted"


_httpd = None



