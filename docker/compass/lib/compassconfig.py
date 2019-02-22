from ConfigParser import *
import pprint
import os
import urllib
import socket
import sys
import subprocess
import shlex

class CompassConf:
    class ToolObj(dict):
        def __init__(self, *args):
            dict.__init__(self, *args)

        def execute(self, source='bin', prepend="", append="", file="", **kwargs):
            '''source must be a config key (eg: path, bin)
               returns stdout,stderr,errcode'''
            cmd = prepend + " " + self[source] + file + " " + append
            logerror.LE.debug("Running command: [{0}]".format(cmd))
            cmdsplit = shlex.split(cmd) #Split command with parameters into a list
            p = subprocess.Popen(
                cmdsplit, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)
            stdout, stderr = p.communicate()
            errcode = p.wait()
            return cmd, stdout, stderr, errcode

        def popen(self, source='bin', prepend="", append="", file="", **kwargs):
            '''source must be a config key (eg: path, bin)
               returns popen object'''
            cmd = prepend + " " + self[source] + file + " " + append
            logerror.LE.debug("Running command: [{0}]".format(cmd))
            cmdsplit = shlex.split(cmd)
            p = subprocess.Popen(cmdsplit, **kwargs)
            p.cmd = cmd
            return p

    def __init__(self, path='/etc/compass.cfg'):
        self.rawconfig = ConfigParser()
        self.rawconfig.read(path)
        self.config = {}
        self.parse()
        self.fixdirs(self.config)
        if "tools" in self.config:
            for i, j in self.config["tools"].items():
                if type(j) == dict:
                    self.config["tools"][i] = CompassConf.ToolObj(j)

    def parse(self):
        sections = self.rawconfig.sections()
        for i in sections:
            if not "." in i:
                self.config[i] = dict(self.rawconfig.items(i))
            else:
                si = i.split(".")
                cfg = self.config
                while si:
                    cfg = cfg.setdefault(si.pop(0), {})
                cfg.update(dict(self.rawconfig.items(i)))

    def __getitem__(self, key):
        return self.config[key]

    def fixdirs(self, dic):
        for i in ['bin', 'path']:
            if i in dic:
                if os.path.exists(dic[i]) and os.path.isdir(dic[i]) and dic[i][-1] != '/':
                    dic[i] += '/'

        for i, j in dic.items():
            if type(j) == dict:
                self.fixdirs(j)

    def check(self, paths=True, hosts=False):
        self.check2([], self.config, paths, hosts)

    def check2(self, sections, dic, paths, hosts):
        if paths:
            for i in ["path", "bin"]:
                if i in dic:
                    if not os.path.exists(dic[i]):
                        logerror.LE.warning(
                            "{0} Path does not exist ({1})".format(sections + [i], dic[i]))
        if hosts:
            if "url" in dic:
                try:
                    urllib.urlopen(dic["url"]).read()
                except:
                    logerror.LE.warning("{0} could not connect to ({1})".format(
                        sections + ['url'], dic['url']))
            if 'host' in dic and 'port' in dic:
                try:
                    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                    s.connect((dic["host"], int(dic["port"])))
                    s.close()
                except:
                    logerror.LE.warning("{0} could not connect to ({1})".format(sections + ['host,port'],
                                                                                dic['host'] + ":" + dic['port']))

        for i, j in dic.items():
            if j.__class__.__name__ in ["dict", "ToolObj"]:
                self.check2(sections + [i], j, paths, hosts)

    def __str__(self):
        return pprint.pformat(self.config)


# Initiate Compass Conf object

configfile = os.environ['COMPASSCFG']
COMPASSCFG = CompassConf(configfile)

import logerror

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print "Usage: {0} [-check] [-print]\n - check = checks paths and hosts\n - print = prints configuration\n".format(
            sys.argv[0])
    else:
        if "-print" in sys.argv:
            print COMPASSCFG
        if "-check" in sys.argv:
            print "Checking paths and hosts"
            COMPASSCFG.check(True, True)
            print "If no Warning Messages everything went OK"
