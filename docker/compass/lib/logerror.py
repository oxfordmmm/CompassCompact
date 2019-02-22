from StringIO import StringIO
import atexit
import datetime
from email.mime.text import MIMEText
import logging
import os
import smtplib
import subprocess
import sys
import traceback

import compassconfig
# from submitlog import *


class CompassLogInfo:
    infoSet = False
    sampleId = ""
    workflowName = ""
    workflowVersion = ""
    processName = ""
    processVersion = ""


# used to identify if the current running process is the logsubmission
# process (so that it can avoid an infinite loop by not logging it's own
# processes)
# logSubmissionScriptName = "submitlog.py"
_compassLogInfo = CompassLogInfo()


# def submitLog(status, comment="", extra=""):
#     logRecord = {'specimenId': _compassLogInfo.sampleId,
#                  'accessionNumber': "",
#                  'processName': _compassLogInfo.processName,
#                  'pipelineName': compassconfig.COMPASSCFG["log"]["pipelinename"],
#                  'pipelineVersion': compassconfig.COMPASSCFG["log"]["pipelineversion"],
#                  'processVersion': _compassLogInfo.processVersion,
#                  'processStatus': status,
#                  'workflowName': _compassLogInfo.workflowName,
#                  'workflowVersion': _compassLogInfo.workflowVersion,
#                  'comment': comment,
#                  'extra': extra}
#     result = writeLogRecord(logRecord)


class SplitStreamHandler(logging.Handler):
    def __init__(self, out=sys.stdout, errout=sys.stderr):
        self.errout = errout
        self.stdout = out
        logging.Handler.__init__(self)

    def emit(self, record):
        if record.msg.__class__.__name__ in ["file", "StringIO"]:
            preffix = ""
            if record.args:
                preffix = record.args[0] + ": "
            rl = logging.getLogger("log.error")
            lno = record.levelno
            for i in record.msg:
                i = i.strip()
                rl._log(lno, preffix + i, [])
        else:
            if record.msg.lower().find("compasslog") == 0:
                if record.msg.find(':') > -1:
                    message = record.msg[10:record.msg.find(':')]
                    comment = record.msg[record.msg.find(':') + 1:]
                else:
                    message = record.msg[10:]
                    comment = ""
                submitLog(message, comment)
            else:
                self.emit2(record)

    def emit2(self, record):
        try:
            msg = self.format(record)

            if record.levelno <= logging.WARNING:
                stream = self.stdout
            # stream=open("/tmp/STDOUT","a")
            else:
                stream = self.errout
            # stream=open("/tmp/STDERR","a")
            fs = "%s\n"

            try:
                if (isinstance(msg, unicode) and
                        getattr(stream, 'encoding', None)):
                    ufs = fs.decode(stream.encoding)
                    try:
                        stream.write(ufs % msg)
                    except UnicodeEncodeError:
                        stream.write((ufs % msg).encode(stream.encoding))
                else:
                    stream.write(fs % msg)
            except UnicodeError:
                stream.write(fs % msg.encode("UTF-8"))

            stream.flush()
        # stream.close() #### CERRAR
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)


def dump_exc():
    """
    Print the usual traceback information, followed by a listing of all the
    local variables in each frame.
    """
    if _compassLogInfo.infoSet and not sys.argv[0][-len(logSubmissionScriptName):] == logSubmissionScriptName:
        logging.info("CompassLogError:" + str(sys.exc_info()[1]))

    log_lines = []
    log_lines.append("============= EXCEPTION RAISED [{0}] ==============".format(
        datetime.datetime.now()))
    log_lines.append("CWD: [{1}] CMD: [{0}]".format(
        " ".join(sys.argv), os.getcwd()))

    _tb = sys.exc_info()[2]
    if not _tb:
        sys.exit(-1)
    while 1:
        if not _tb.tb_next:
            break
        _tb = _tb.tb_next
    stack = []
    f = _tb.tb_frame
    while f:
        stack.append(f)
        f = f.f_back
    stack.reverse()
    for i in traceback.format_exc().strip().split("\n"):
        log_lines.append("EXCEPTION: " + i)

    # logging.critical( "EXCEPTION: "+traceback.format_exc().strip().replace("\n","\nEXCEPTION: "))

    for frame in stack:
        log_lines.append(
            "\tVARS: Frame %s in %s at line %s" % (frame.f_code.co_name, frame.f_code.co_filename, frame.f_lineno))
        for key, value in frame.f_locals.items():
            if key == "self" or key.startswith("__"):
                continue
            if sys.getsizeof(str(value)) > 500:
                value = "[Too large value]"

            try:
                log_lines.append("\tVARS: |_\t{0} = {1}".format(key, value))
            except:
                pass
    log_lines.append("============= END EXCEPTION [{0}] ==============".format(
        datetime.datetime.now()))

    # Print all the lines in the log_lines to critical log.
    for line in log_lines:
        logging.critical(line)

        # Email all the lines in log_lines to an email.
    #     notify_error("\r\n".join(log_lines))
    # with open("/Users/thanhlv/Pipeline/compass/experiments/log.txt", "w") as logfile:
    #     logfile.write(log_lines)
        
    sys.exit(-1)


def notify_error(message, email="mmmsequencing@ndm.ox.ac.uk"):
    """Email the error to an account."""

    # Create a text/plain message
    msg = MIMEText(message)

    msg['Subject'] = "ERROR"
    msg['To'] = email
    # FIXME: Don't know if From required.
    #     msg["From"] = email

    # envelope header.
    s = smtplib.SMTP('localhost')
    s.sendmail(email, [email], msg.as_string())
    s.quit()


# look for compass log-specific data in the argument list, and remove when found
# structure allows for process name and version to be ommitted, but not
# workflow name and version
def setCompassLogInfo():
    if '-wfn' in sys.argv and '-wfv' in sys.argv:
        if '-u' in sys.argv:
            idIndex = sys.argv.index('-u')
            _compassLogInfo.sampleId = sys.argv[idIndex + 1]
            sampleIdProvided = True
        elif '-ux' in sys.argv:
            idIndex = sys.argv.index('-ux')
            sys.argv.pop(idIndex)
            _compassLogInfo.sampleId = sys.argv.pop(idIndex)
            sampleIdProvided = True
        else:
            sampleIdProvided = False
        if sampleIdProvided:
            nameindex = sys.argv.index('-wfn')
            sys.argv.pop(nameindex)
            _compassLogInfo.workflowName = sys.argv.pop(nameindex)
            versionindex = sys.argv.index('-wfv')
            sys.argv.pop(versionindex)
            _compassLogInfo.workflowVersion = sys.argv.pop(versionindex)
            if '-prn' in sys.argv and '-prv' in sys.argv:
                processNameIndex = sys.argv.index('-prn')
                sys.argv.pop(processNameIndex)
                _compassLogInfo.processName = sys.argv.pop(processNameIndex)
                processVersionIndex = sys.argv.index('-prv')
                sys.argv.pop(processVersionIndex)
                _compassLogInfo.processVersion = sys.argv.pop(
                    processVersionIndex)
            _compassLogInfo.infoSet = True


setCompassLogInfo()

LE = logging.getLogger('log.error')
LE.addHandler(SplitStreamHandler())
LE.setLevel(logging.DEBUG)
LE.handlers[0].setFormatter(
    logging.Formatter("** %(levelname)8s:%(filename)s:%(funcName)s[%(lineno)d]:%(asctime)s :: %(message)s"))

LE.info("CMD: [ {0} ]".format(" ".join(sys.argv)))
LE.info("CWD: [ {0} ]".format(os.getcwd()))

if _compassLogInfo.infoSet and not sys.argv[0][-len(logSubmissionScriptName):] == logSubmissionScriptName:
    LE.info("CompassLogStart")


def endPythonScript():
    LE.info("END-CMD: [ {0} ]".format(" ".join(sys.argv)))
    if _compassLogInfo.infoSet and not sys.argv[0][-len(logSubmissionScriptName):] == logSubmissionScriptName:
        LE.info("CompassLogFinish")


atexit.register(endPythonScript)

if __name__ == "__main__":
    LE.debug('debug message (this is stdout)')
    LE.info('info message (this is stdout)')
    LE.warn('warn message (this is stdout)')
    LE.error('error message (this is stderr)')
    LE.critical('critical message (this is stderr)')

    import code

    code.interact(local=locals())

    class a:
        def st(self):
            self.qt()

        def qt(self):
            qwe = 2
            ewr = er

    b = a()
    try:
        b.st()
    except:
        dump_exc()
