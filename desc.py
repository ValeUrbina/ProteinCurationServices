from cherrypy.lib.static import serve_file
import glob
import os.path

import cherrypy
from cherrypy.lib.static import serve_file

# ver como pasar el archivo al front: AUXILIO


class Root:
    def index(self, directory="."):
        html = """<html><body><h2>Here are the files in the selected directory:</h2>
        <a href="index?directory=%s">Up</a><br />
        """ % os.path.dirname(os.path.abspath(directory))

        for filename in glob.glob(directory + '/*'):
            absPath = os.path.abspath(filename)
            if os.path.isdir(absPath):
                html += '<a href="/index?directory=' + absPath + \
                    '">' + os.path.basename(filename) + "</a> <br />"
            else:
                html += '<a href="/download/?filepath=' + absPath + \
                    '">' + os.path.basename(filename) + "</a> <br />"

        html += """</body></html>"""
        return html
    index.exposed = True


class Download:
    def index(self, filepath):
        index.exposed = True
        return serve_file(filepath, "application/x-download", "attachment")


if __name__ == '__main__':
    root = Root()
    root.download = Download()
    cherrypy.quickstart(root)
