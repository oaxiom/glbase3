
"""

a simple progress bar for some of the longer functions expected to take a long-time.

Inspired from the progress bar in urlgrabber.

"""



import sys

from . import config

class progressbar:
    def __init__(self, maximum_value, output=sys.stderr):
        """
        Initialise the progress bar

        **Arguments**

            maximum_value (Required)
                the maximum_value to move the bar up to.

            output (Optional, defaults to stderr)
                the output device to use.
        """
        self.maximum = maximum_value
        if maximum_value <= 0: # special case, probably being sent a delayedlist
            self.maximum = -1

        self.__writer = output
        self.__barwidth = 30 # bar_width in characters. This may need to change on later computers
                             # with larger terminals
        self.__last_percent = -1 # only print if __last_pecent is incremented.

    def update(self, new_value):
        """
        Update progress meter with new_value

        **Arguments**

            new_value (Required)
                should be some number between 0 .. maximum_value
        """
        if self.maximum == -1:
            return(False) # disable progress bar in wierd situations

        t_percent_done = int(((new_value+1) / self.maximum) * self.__barwidth)

        if t_percent_done > self.__last_percent:
            percent_done = int(((new_value+1) / self.maximum) *100)

            bar = "".join(["=" for x in range(t_percent_done)] + ["-" for x in range(self.__barwidth-t_percent_done)])
            if not config.SILENT: self.__writer.write("\r[%s] %s%% (%s/%s)" % (bar, percent_done, new_value+1, self.maximum))
            self.__last_percent = t_percent_done

        if new_value+1 >= self.maximum: # if the last line, reset the console so the result overprints the progress bar.
            if not config.SILENT: self.__writer.write("\r") # pad out to overwrite the previous bar.
            if not config.SILENT: self.__writer.write("\r                                                        ") # pad out to overwrite the previous bar.
            if not config.SILENT: self.__writer.write("\r") # pad out to overwrite the previous bar.

if __name__ == "__main__":
    """
    A quick tester.
    """
    import time

    p = progressbar(100)
    for i in range(100):
        time.sleep(0.1)
        p.update(i)
