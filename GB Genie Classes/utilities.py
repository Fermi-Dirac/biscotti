import sys
import os


def tabulate_item(row, col_widths, sep="  "):
    """Convert a string list representing a row into a formatted string.

    Args:
        row (str list): A list of strings representing a row in a table, each 
            item represents a column within a row of the table.
        col_widths (int list): A list of integers representing the column 
            width of the table.
        sep (str, optional): The separation character between rows. Default is
            two spaces.

    Returns:
        str: A string representing a tabulated row.
    """
    str_reps = []
    row_wid = len(row)
    for i in range(0, len(row)):
        wid = col_widths[i]
        str_reps.append(row[i].ljust(wid))
    out = sep.join(str_reps)
    return out


def tabulate(rows, sep="  "):
    """Generate a formatted string representing a table.

    Args:
        rows (str list list): A list in which each element is a list of strings
            representing a row in the table.
        sep (str, optional): The separation character between rows. Default is
            two spaces.

    Returns:
        str: Formatted string representing a table.
    """
    cols = max(map(len, rows))
    col_widths = []
    for i in range(0, cols):
        w = 0
        for c in rows:
            if i < len(c):
                w = max(w, len(c[i]))
        col_widths.append(w)
    res = "\n".join(map(lambda x: tabulate_item(x, col_widths, sep), rows))
    return res


def open_read_file(path, extension):
    if (path.split('.')[-1] != extension):
        raise ValueError('File %s is not of extension %s.' % (path, extension))
    if (not os.path.isfile(path)):
        raise ValueError('File %s does not exist.' % path)
    return open(path, 'r')


def open_write_file(path, overwrite_protect=True):
    if False: #(os.path.isfile(path) and overwrite_protect):
        print('File %s exists.' % path)
        decision = raw_input('Type "Y" if you want to overwrite or' +
                             ' any other character(s) to quit.\n')
        if (decision != 'Y'):
            sys.exit('You have aborted the task.')
    return open(path, 'w')
