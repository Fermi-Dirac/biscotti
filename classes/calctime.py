import re
import datetime as dt
from collections import OrderedDict as odict
import logging
# Setup a logger
logger = logging.getLogger(__name__)
loglevel = logging.INFO
logger.setLevel(loglevel)
console_handler = logging.StreamHandler()
console_handler.setLevel(loglevel)
console_handler.setFormatter(logging.Formatter('%(name)s - %(levelname)s - %(message)s'))
logger.addHandler(console_handler)
# End logging.

try:
    import matplotlib.pyplot as plt
    has_mpl = True
    plt.style.use('fivethirtyeight')
except ImportError:
    logger.error("Cannot load matplotlib! Plotting disabled")
    has_mpl = False

class CalcTime(object):
    def __init__(self, startdt = None, enddt = None, timedict=None, subprocdict=None, procorder = None):
        if timedict is None:
            timedict = {}
        if subprocdict is None:
            subprocdict = {}
        if procorder is None:
            procorder = []
        self.startdt = startdt
        self.enddt = enddt
        self.timedict = timedict
        self.subprocdict = subprocdict
        self.procorder = procorder

    @staticmethod
    def get_time(outputfile):
        logger.debug("Starting calctime.get_time on " + outputfile)
        regexstart = r'[ ]*init_run[\s\:]+[0-9\.s]+'
        regexend = r'This run was terminated on:[\s]+'
        startson_dtformat = '%d%b%Y at %H:%M:%S'
        starttime_regex  = r'[0-9]+[A-Za-z]{3}[0-9]{4} at[\s]+[0-9\s]+:[0-9\s]+:[\s]*[0-9]+'
        ends_dtformat = '%H:%M:%S  %d%b%Y'
        endtime_regex = r'[0-9]+:[0-9\s]+:[0-9\s]+[A-Za-z]{3}[0-9]{4}'

        with open(outputfile, 'r') as fileobj:

            filestring = fileobj.read()
            startresult = re.search(regexstart, filestring)
            endresult = re.search(regexend, filestring)

            if startresult and endresult:
                logger.debug("Found start at " + str(startresult.start()))
                calctimestring = filestring[startresult.start(): endresult.end()]
            else:
                calctimestring = ''
            logger.debug("Match follows\n" + calctimestring)
            rows = calctimestring.split('\n')
            timedict = odict()
            subprockey = 'Main'
            subprocdict = {subprockey: []}
            procorder = [subprockey]
            for row in rows:
                logger.debug("Now on row: " + row)
                if row.strip() is '':
                    subprockey = None
                else:
                    keyresult = re.search(r'Called by [A-Za-z\*_]+', row)
                    keyresult2 = re.search(r'[A-Za-z]+\sroutines', row)
                    if keyresult or keyresult2:
                        logger.debug("New key to be " + row)
                        if keyresult:
                            subprockey = keyresult.group(0).split('by')[1].strip()
                        else:
                            subprockey = keyresult2.group(0).strip()
                        procorder.append(subprockey)
                        subprocdict[subprockey] = []
                        logger.debug("New subprocess key")
                    else:
                        if subprockey is None:
                            pass
                            # Here we could add code to parse the PWSCF step at the end giving total time
                        else:
                            procname = row[5:18].strip()
                            values = [float(row[19:29]), float(row[34:44]), float(row[52:60])]
                            # for i, value in enumerate(values):
                            #     logger.debug(value + " : converting to float")
                            #     values[i] = float(re.compile(r'[^\d.]+').sub('', value))  #strip off that annoying 's'
                            timedict[procname] = values
                            subprocdict[subprockey].append(procname)
                            logger.debug("New entry of process: " + procname + " with values: " + str(values))
            # Now get the start and end datetimes
            start_dt_result = re.search(starttime_regex, filestring)
            if start_dt_result:
                dtsplit = start_dt_result.group(0).split(' at ')
                dtstring = dtsplit[0] + ' at ' + dtsplit[1].replace(' ', '0')
                startdt = dt.datetime.strptime(dtstring, startson_dtformat)
            else:
                startdt = dt.datetime.min

            end_dt_result = re.search(endtime_regex, filestring)
            if end_dt_result:
                dtsplit = end_dt_result.group(0).split('  ')
                dtstring = dtsplit[0].replace(' ', '0') + '  ' + dtsplit[1]
                enddt = dt.datetime.strptime(dtstring, ends_dtformat)
            else:
                enddt = dt.datetime.max

        return CalcTime(startdt, enddt, timedict, subprocdict, procorder)

    def pie_charts(self, type = 'CPU'):
        if not has_mpl:
            logger.error("Error, matplotlib not loaded!")
            return None
        timeindex = {'CPU': 0 , 'cpu': 0, 'WALL' : 1, 'wall' : 1, 'calls' : 2}
        if type not in timeindex:
            type = 'CPU'
        title = type + ' time breakdown '
        fig = plt.figure()
        fig.suptitle(title, fontsize = 14, fontweight = 'bold')
        for i, subproc in enumerate(self.procorder):
            labels = []
            sizes = []
            for key in self.subprocdict[subproc]:
                units = "(" + "%.2f" % (self.timedict[key][timeindex[type]] / (60 * 60)) + " Hrs)"
                if type == 'calls':
                    units = "(" + "%.2f" % (self.timedict[key][timeindex[type]]) + " calls)"
                labels.append(key + " " + units)
                sizes.append(self.timedict[key][0])
            if len(sizes) > 0:
                plt.subplot(2, len(self.procorder) / 2 + 1, i + 1)
                # fig = plt.figure(i+1)
                # ax = fig.gca()
                siz_lab = sorted(zip(sizes, labels), key= lambda tup: tup[0], reverse = True)
                logger.debug("New pie charge, sorted zip is: " + str(siz_lab))
                pie_size = sum(sizes)
                new_sizes, new_labels = map(list, zip(*[tup for tup in siz_lab if tup[0]/pie_size > 0.01])) # Unzip
                if len(new_sizes) < len(sizes): # if we trimmed some
                    other_size, other_labels = map(list, zip(*[tup for tup in siz_lab if tup[0] / pie_size < 0.01]))
                    new_sizes.append(sum(other_size)/pie_size)
                    new_labels.append('Others: ' + '\n'.join(other_labels))
                plt.pie(new_sizes, labels=new_labels, shadow=True, startangle=45, autopct='%1.1f %%')
                plt.axis('equal')
                plt.title(subproc)
        return fig

    def pie_chart(self, axes=None, type='cpu', subproc = None):
        if has_mpl is False:
            logger.error("Matplotlib dependency not loaded. Plotting disabled!")
            return None
        if axes is None:
            axes = plt.gca()
        timeindex = {'CPU': 0, 'cpu': 0, 'WALL': 1, 'wall': 1, 'calls': 2}
        if type not in timeindex:
            type = 'CPU'
        if subproc is None:
            subproc = self.procorder[0]
        title = type + ' time spent on ' + str(subproc)
        labels = []
        sizes = []
        for key in self.subprocdict[subproc]:
            if type == 'calls':
                units = "(" + "%.2f" % (self.timedict[key][timeindex[type]]) + " calls)"
            else:
                units = "(" + "%.2f" % (self.timedict[key][timeindex[type]] / (60 * 60)) + " Hrs)"
            labels.append(key + " " + units)
            sizes.append(self.timedict[key][0])
        if len(sizes) > 0:
            siz_lab = sorted(zip(sizes, labels), key=lambda tup: tup[0], reverse=True)
            logger.debug("New pie charge, sorted zip is: " + str(siz_lab))
            pie_size = sum(sizes)
            new_sizes, new_labels = map(list, zip(*[tup for tup in siz_lab if tup[0] / pie_size > 0.01]))  # Unzip
            if len(new_sizes) < len(sizes):  # if we trimmed some
                other_size, other_labels = map(list, zip(*[tup for tup in siz_lab if tup[0] / pie_size < 0.01]))
                new_sizes.append(sum(other_size) / pie_size)
                new_labels.append('Others: ' + '\n'.join(other_labels))
            axes.pie(new_sizes, labels=new_labels, shadow=True, startangle=45, autopct='%1.1f %%')
            axes.axis('equal')
            axes.set_title(title)
        return axes