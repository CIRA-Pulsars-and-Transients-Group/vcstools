import sys

def progress_bar(it, prefix="", size=60, file=sys.stdout):
    """
    I stole this code from here: https://stackoverflow.com/questions/3160699/python-progress-bar
    It's used like this: for i in progressbar(range(15), "Computing: ", 40):
    """
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