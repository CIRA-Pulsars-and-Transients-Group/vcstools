import sys

def progress_bar(it, prefix="", size=60, file=sys.stdout):
    """I stole this code from here: https://stackoverflow.com/questions/3160699/python-progress-bar

    Parameters
    ----------
    it : `list`
        The list to iterate over.
    prefix : `str`
        The prefix do display in the progress bar. |br| Default: "".
    size : `int`
        The length of the progress bar to display in characters. |br| Default: 60.
    file : stdout
        The output of the progress bar. |br| Default: `sys.stdout`.
    
    Examples
    --------
    >>> for i in progressbar(range(15), "Computing: ", 40):
    """
    it = list(it)
    count = len(it)
    def show(j):
        x = int(size*j/count)
        file.write("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x), j, count))
        file.flush()
    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    file.write("\n")
    file.flush()